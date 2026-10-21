# hash:sha256:67309b0a90db79eaa06a3936866ba7603a087477b48919808cbd7625c6dc674a
FROM registry.codeocean.ir-scc-codeocean-prod.awscloud.abbvienet.com/codeocean/cpu-qmg-oncology-co-red-repo-base:v1.0.0

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
