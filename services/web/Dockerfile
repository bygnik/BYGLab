# Multi-stage: build the wasm solver, build React, serve with nginx.
# Preserves the repo's real directory nesting (rust/ as a sibling of
# services/web/) so the relative import in src/wasm-worker.js resolves
# identically here and in local dev.

FROM rust:1-slim AS wasm-build
RUN rustup target add wasm32-unknown-unknown && cargo install wasm-pack
WORKDIR /repo/rust
COPY rust/ ./
RUN wasm-pack build crates/byglab-wasm --target web --release

FROM node:20-alpine AS build
WORKDIR /repo/services/web
COPY services/web/package.json services/web/package-lock.json* ./
RUN npm install
COPY services/web/ ./
COPY --from=wasm-build /repo/rust/crates/byglab-wasm/pkg /repo/rust/crates/byglab-wasm/pkg
RUN npm run build

FROM nginx:1.27-alpine
COPY services/web/nginx.conf /etc/nginx/conf.d/default.conf
COPY --from=build /repo/services/web/dist /usr/share/nginx/html
EXPOSE 80
