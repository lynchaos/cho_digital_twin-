---
name: Express + http-proxy-middleware POST body fix
description: http-proxy-middleware v4 fails to forward POST bodies if express.json() runs first
---

# Rule
Do NOT place `express.json()` (or any body-parsing middleware) before `createProxyMiddleware()` in Express.

**Why:** Body parsers consume the Node.js request stream. Once read, the stream is empty. http-proxy-middleware v4 tries to pipe the raw stream to the upstream target — if it's already been consumed, the upstream gets no body. The proxy hangs until timeout. GET requests work fine because they have no body.

**How to apply:**
- Place the proxy middleware BEFORE `express.json()` / `express.urlencoded()`.
- The `/api` router (which needs parsed JSON) mounts AFTER the proxy block.
- This ordering lets the proxy forward raw stream bytes for proxied routes while body parsing still works for direct Express routes.

**Symptom:** POST through proxy times out ("request aborted" after N seconds) even though direct calls to the upstream port succeed instantly. GETs and POLLs work fine.
