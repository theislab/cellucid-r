// The front end of the Cellucid notebook widget, an anywidget module shared
// byte for byte by cellucid-python (cellucid/jupyter/static/widget.js, behind
// cellucid.Widget) and cellucid-r (inst/widget/widget.js, behind
// cellucid_view() in a Jupyter notebook). It frames the Cellucid app at
// `${web_source_url}/?viewer=embedded` and relays between the page and the
// kernel:
//
// - The page says `ready` to its parent. Only a message from this output's
//   own frame and the app's origin is accepted; the answer, `connect`, is
//   posted to that exact origin and carries one MessagePort and the viewer's
//   id and token, which this module asks the kernel for and keeps in memory
//   only (never in a trait, an address, storage or a log).
// - Everything after travels on that port: the viewer channel's frames
//   (`open`, `text`, `binary`, `close`, by socket `s`) and the page's
//   requests (`fetch`), each relayed to the kernel as a widget message named
//   by this output (`v`), and the kernel's frames and answers back.
// - Requests that arrive together go to the kernel as one message; a binary
//   frame goes in pieces, so no message is larger than a notebook carries.
// - The kernel answers its messages one after another, and none while a
//   cell runs: at most two batches of requests are with it at once (one
//   being answered, the next ready), and the rest wait here, where a request
//   the page withdraws (`cancel`) is dropped before the kernel ever reads
//   for it. A withdrawn request already with the kernel is answered, and its
//   answer dropped by the page.
// - A kernel whose comms carry no binary buffers (IRkernel) sets the trait
//   `_binary` to "base64": binary bodies then travel as base64 text (`b64`).
//
// The page's side is cellucid/assets/js/data/viewer-embed.js in the web app.

const EMBED_KEY = 'cellucidViewerEmbed';

// A binary frame is sent to the kernel in pieces of this size (the
// notebook server's connection refuses very large messages).
const PIECE_BYTES = 1024 * 1024;

// How long after the app's page has loaded the output waits for its `ready`
// before it says what may be wrong (it keeps waiting).
const SILENT_PAGE_NOTICE_MS = 15000;

function randomId() {
  const bytes = new Uint8Array(16);
  crypto.getRandomValues(bytes);
  return Array.from(bytes, byte => byte.toString(16).padStart(2, '0')).join('');
}

/** An ArrayBuffer of its own holding a buffer the kernel sent. */
function ownBytes(view) {
  if (view instanceof ArrayBuffer) return view.slice(0);
  return view.buffer.slice(view.byteOffset, view.byteOffset + view.byteLength);
}

function fromBase64(text) {
  const binary = atob(text);
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index++) bytes[index] = binary.charCodeAt(index);
  return bytes.buffer;
}

function toBase64(bytes) {
  let binary = '';
  for (let start = 0; start < bytes.length; start += 0x8000) {
    binary += String.fromCharCode.apply(null, bytes.subarray(start, start + 0x8000));
  }
  return btoa(binary);
}

function render({ model, el }) {
  const view = randomId();
  const appUrl = model.get('web_source_url');
  const base64 = model.get('_binary') === 'base64';
  const appOrigin = new URL(appUrl).origin;

  const box = document.createElement('div');
  box.className = 'cellucid-widget';
  box.style.height = `${model.get('height')}px`;
  const frame = document.createElement('iframe');
  frame.className = 'cellucid-widget-frame';
  frame.title = 'Cellucid viewer';
  frame.allow = 'clipboard-write; fullscreen';
  frame.referrerPolicy = 'no-referrer';
  const status = document.createElement('div');
  status.className = 'cellucid-widget-status';
  status.setAttribute('role', 'status');
  status.textContent = `Loading Cellucid from ${appOrigin}…`;
  box.append(frame, status);
  el.append(box);

  let port = null;
  let credentials = null;
  let pageReady = false;
  let ended = false;
  let silentNotice = null;
  /** Requests waiting to go to the kernel together. */
  let requests = [];
  let flushing = false;
  /** The batches with the kernel: the ids each still waits an answer for. */
  const outstanding = [];

  const say = text => {
    status.textContent = text;
    status.hidden = text === null;
  };

  const sendKernel = (content, buffers) => {
    if (!ended) model.send({ ...content, v: view }, undefined, buffers);
  };

  const flushRequests = () => {
    flushing = false;
    if (requests.length === 0 || outstanding.length >= 2) return;
    sendKernel({ t: 'fetch', requests });
    outstanding.push(new Set(requests.map(request => request.id)));
    requests = [];
  };

  /** The kernel answered `id`: a batch it completes lets the next one go. */
  const answered = id => {
    const batch = outstanding.findIndex(ids => ids.delete(id));
    if (batch >= 0 && outstanding[batch].size === 0) {
      outstanding.splice(batch, 1);
      flushRequests();
    }
  };

  const fromPage = message => {
    if (message === null || typeof message !== 'object') return;
    switch (message.t) {
      case 'fetch':
        requests.push({ id: message.id, method: message.method, path: message.path, range: message.range ?? null });
        if (!flushing) {
          flushing = true;
          setTimeout(flushRequests, 0);
        }
        return;
      case 'cancel':
        // Not with the kernel yet: never sent.
        requests = requests.filter(request => request.id !== message.id);
        return;
      case 'open':
        sendKernel({ t: 'open', s: message.s });
        return;
      case 'text':
        sendKernel({ t: 'text', s: message.s, text: message.text });
        return;
      case 'binary': {
        const bytes = new Uint8Array(message.bytes);
        const sendPiece = (piece, last) => {
          if (base64) sendKernel({ t: 'binary', s: message.s, last, b64: toBase64(piece) });
          else sendKernel({ t: 'binary', s: message.s, last }, [piece.slice().buffer]);
        };
        if (bytes.byteLength === 0) {
          sendPiece(bytes, true);
          return;
        }
        for (let start = 0; start < bytes.byteLength; start += PIECE_BYTES) {
          const end = Math.min(start + PIECE_BYTES, bytes.byteLength);
          sendPiece(bytes.subarray(start, end), end === bytes.byteLength);
        }
        return;
      }
      case 'close':
        sendKernel({ t: 'close', s: message.s, code: message.code });
        return;
      default:
    }
  };

  const connect = () => {
    if (!pageReady || credentials === null || port !== null || ended) return;
    const channel = new MessageChannel();
    port = channel.port1;
    port.onmessage = event => fromPage(event.data);
    frame.contentWindow.postMessage(
      { [EMBED_KEY]: 'connect', viewerId: credentials.viewerId, viewerToken: credentials.viewerToken },
      appOrigin,
      [channel.port2]
    );
    say(null);
  };

  const onWindowMessage = event => {
    if (event.source !== frame.contentWindow || event.origin !== appOrigin) return;
    if (event.data === null || typeof event.data !== 'object' || event.data[EMBED_KEY] !== 'ready') return;
    clearTimeout(silentNotice);
    if (port !== null) {
      // The page loaded again: what it had open is gone, and it gets a new port.
      sendKernel({ t: 'detach' });
      port.close();
      port = null;
      // What the page had asked went with it.
      requests = [];
      outstanding.length = 0;
    }
    pageReady = true;
    if (credentials === null) say('Waiting for the notebook kernel…');
    connect();
  };

  const fromKernel = (message, buffers) => {
    if (message === null || typeof message !== 'object' || message.v !== view) return;
    switch (message.t) {
      case 'credentials':
        credentials = { viewerId: message.viewerId, viewerToken: message.viewerToken };
        connect();
        return;
      case 'answers':
        for (const answer of message.answers) {
          const bytes = typeof answer.b64 === 'string'
            ? fromBase64(answer.b64)
            : answer.buffer >= 0 ? ownBytes(buffers[answer.buffer]) : new ArrayBuffer(0);
          port?.postMessage(
            { t: 'answer', id: answer.id, status: answer.status, headers: answer.headers, bytes },
            [bytes]
          );
          answered(answer.id);
        }
        return;
      case 'text':
        port?.postMessage({ t: 'text', s: message.s, text: message.text });
        return;
      case 'binary': {
        const bytes = typeof message.b64 === 'string' ? fromBase64(message.b64) : ownBytes(buffers[0]);
        port?.postMessage({ t: 'binary', s: message.s, bytes }, [bytes]);
        return;
      }
      case 'close':
        port?.postMessage({ t: 'close', s: message.s, code: message.code, reason: message.reason ?? '' });
        return;
      default:
    }
  };

  const onHeight = () => {
    box.style.height = `${model.get('height')}px`;
  };

  frame.addEventListener('load', () => {
    if (pageReady) return;
    clearTimeout(silentNotice);
    silentNotice = setTimeout(() => {
      if (pageReady || ended) return;
      say(
        `The page at ${appUrl} has not answered as a Cellucid viewer. The notebook's ` +
        'browser loads the app from there, so it needs to reach that address (or web_source_url ' +
        'set to a copy you host); an older build of the app cannot be embedded either.'
      );
    }, SILENT_PAGE_NOTICE_MS);
  });
  window.addEventListener('message', onWindowMessage);
  model.on('msg:custom', fromKernel);
  model.on('change:height', onHeight);
  frame.src = `${appUrl}/?viewer=embedded`;
  sendKernel({ t: 'attach' });

  return () => {
    sendKernel({ t: 'detach' });
    ended = true;
    clearTimeout(silentNotice);
    window.removeEventListener('message', onWindowMessage);
    model.off('msg:custom', fromKernel);
    model.off('change:height', onHeight);
    port?.close();
    port = null;
  };
}

export default { render };
