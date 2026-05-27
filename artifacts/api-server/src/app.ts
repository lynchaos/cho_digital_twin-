import express, { type Express } from "express";
import cors from "cors";
import pinoHttp from "pino-http";
import { createProxyMiddleware } from "http-proxy-middleware";
import router from "./routes";
import { logger } from "./lib/logger";

const app: Express = express();

app.use(
  pinoHttp({
    logger,
    serializers: {
      req(req) {
        return {
          id: req.id,
          method: req.method,
          url: req.url?.split("?")[0],
        };
      },
      res(res) {
        return {
          statusCode: res.statusCode,
        };
      },
    },
  }),
);
app.use(cors());

// Proxy /gem/* → Python FastAPI GEM Reduction Service on port 8082
// IMPORTANT: do NOT put express.json() / bodyParser before this proxy —
// consuming the request body stream breaks http-proxy-middleware v4 POST forwarding.
app.use(
  createProxyMiddleware({
    target: "http://localhost:8082",
    changeOrigin: true,
    pathFilter: "/gem",
    proxyTimeout: 300_000,
    timeout: 300_000,
    on: {
      error: (err, _req, res) => {
        logger.warn({ err }, "GEM service proxy error");
        const httpRes = res as express.Response;
        if (typeof httpRes.status === "function" && !httpRes.headersSent) {
          httpRes.status(503).json({
            error: "GEM Reduction Service unavailable",
            detail: "The Python backend is starting up or encountered an error. Try again in a few seconds.",
          });
        }
      },
    },
  }),
);

app.use("/api", router);

export default app;
