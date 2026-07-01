# Single Dockerfile: builder stage builds RDKit; final stage is Debian Trixie runtime
ARG PG_MAJOR=18
ARG RDKIT_VERSION=Release_2025_09_1

########################################
# Stage: builder (build RDKit)
########################################
FROM debian:trixie AS builder

ARG PG_MAJOR
ARG RDKIT_VERSION
ENV DEBIAN_FRONTEND=noninteractive
ENV CMAKE_PREFIX_PATH=/usr/local

# Install build deps and PGDG repo so postgresql-18 is available during builder
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates wget gnupg2 lsb-release dirmngr && \
    CODENAME="$(lsb_release -cs)"; \
    wget -qO /usr/share/keyrings/pgdg.gpg https://www.postgresql.org/media/keys/ACCC4CF8.asc; \
    echo "deb [signed-by=/usr/share/keyrings/pgdg.gpg] http://apt.postgresql.org/pub/repos/apt ${CODENAME}-pgdg main" \
      > /etc/apt/sources.list.d/pgdg.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
      build-essential cmake git ca-certificates curl wget python3 python3-pip \
      postgresql-${PG_MAJOR} postgresql-server-dev-${PG_MAJOR} \
      libeigen3-dev libcairo2-dev libfreetype6-dev libpng-dev libjpeg-dev \
      zlib1g-dev pkg-config libboost-all-dev && \
    rm -rf /var/lib/apt/lists/*

# Fetch RDKit and build
WORKDIR /opt
RUN git clone --depth 1 --branch "${RDKIT_VERSION}" https://github.com/rdkit/rdkit.git /opt/rdkit

WORKDIR /opt/rdkit_build
RUN set -eux; \
    PG_INC="$(/usr/lib/postgresql/${PG_MAJOR}/bin/pg_config --includedir)"; \
    PG_INC_SRV="$(/usr/lib/postgresql/${PG_MAJOR}/bin/pg_config --includedir-server)"; \
    PG_LIB="$(/usr/lib/postgresql/${PG_MAJOR}/bin/pg_config --libdir)"; \
    cmake -S /opt/rdkit -B /opt/rdkit_build \
      -DCMAKE_BUILD_TYPE=Release \
      -DRDK_BUILD_PGSQL=ON \
      -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
      -DRDK_BUILD_INCHI_SUPPORT=ON \
      -DRDK_BUILD_AVALON_SUPPORT=ON \
      -DRDK_BUILD_CAIRO_SUPPORT=ON \
      -DRDK_PGSQL_STATIC=ON \
      -DPostgreSQL_CONFIG=/usr/lib/postgresql/${PG_MAJOR}/bin/pg_config \
      -DPG_CONFIG=/usr/lib/postgresql/${PG_MAJOR}/bin/pg_config \
      -DPostgreSQL_INCLUDE_DIR="${PG_INC}" \
      -DPostgreSQL_TYPE_INCLUDE_DIR="${PG_INC_SRV}" \
      -DPostgreSQL_LIBRARY_DIR="${PG_LIB}"; \
    cmake --build /opt/rdkit_build -j"$(nproc)"; \
    cmake --install /opt/rdkit_build --strip

# Preserve artifacts in canonical places for copying into final
RUN mkdir -p /usr/lib/postgresql/${PG_MAJOR}/lib /usr/share/postgresql/${PG_MAJOR}/extension /usr/local/lib && \
    cp -a /usr/lib/postgresql/${PG_MAJOR}/lib/rdkit.so /usr/lib/postgresql/${PG_MAJOR}/lib/ 2>/dev/null || true && \
    cp -a /usr/share/postgresql/${PG_MAJOR}/extension/rdkit* /usr/share/postgresql/${PG_MAJOR}/extension/ 2>/dev/null || true && \
    cp -a /usr/lib/x86_64-linux-gnu/libboost_*.so* /usr/local/lib/ 2>/dev/null || true

########################################
# Stage: final (runtime on Debian Trixie)
########################################
FROM debian:trixie AS final

ARG PG_MAJOR
ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/lib/postgresql/${PG_MAJOR}/bin:${PATH}"
ENV PGDATA=/var/lib/postgresql/data

# Add PGDG repo in final stage so postgresql-18 packages are available
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates wget gnupg2 lsb-release dirmngr && \
    CODENAME="$(lsb_release -cs)"; \
    wget -qO /usr/share/keyrings/pgdg.gpg https://www.postgresql.org/media/keys/ACCC4CF8.asc; \
    echo "deb [signed-by=/usr/share/keyrings/pgdg.gpg] http://apt.postgresql.org/pub/repos/apt ${CODENAME}-pgdg main" \
      > /etc/apt/sources.list.d/pgdg.list && \
    apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates wget gnupg2 lsb-release dirmngr gettext procps \
      libcairo2 libfreetype6 libpng16-16 libjpeg62-turbo zlib1g gosu \
      postgresql-${PG_MAJOR} postgresql-server-dev-${PG_MAJOR} && \
    rm -rf /var/lib/apt/lists/*

# Copy artifacts from the local builder stage (NOT from a deleted external image)
COPY --from=builder /usr/lib/postgresql/${PG_MAJOR}/lib/rdkit.so /usr/lib/postgresql/${PG_MAJOR}/lib/
COPY --from=builder /usr/share/postgresql/${PG_MAJOR}/extension/rdkit* /usr/share/postgresql/${PG_MAJOR}/extension/
COPY --from=builder /usr/local/lib/ /usr/local/lib/

# change UID/GID and fix ownership (optional)
RUN set -eux; \
    OLD_UID="$(id -u postgres)"; OLD_GID="$(id -g postgres)"; \
    usermod -u 999 postgres || true; groupmod -g 999 postgres || true; \
    find / -path /proc -prune -o -user "${OLD_UID}" -exec chown -h 999 {} + 2>/dev/null || true; \
    find / -path /proc -prune -o -group "${OLD_GID}" -exec chgrp -h 999 {} + 2>/dev/null || true

# Make common PG tools available in /usr/bin
RUN ln -sf /usr/lib/postgresql/${PG_MAJOR}/bin/initdb /usr/bin/initdb || true && \
    ln -sf /usr/lib/postgresql/${PG_MAJOR}/bin/pg_isready /usr/bin/pg_isready || true && \
    ln -sf /usr/lib/postgresql/${PG_MAJOR}/bin/pg_ctl /usr/bin/pg_ctl || true && \
    ln -sf /usr/lib/postgresql/${PG_MAJOR}/bin/psql /usr/bin/psql || true

# Vendor entrypoint (place official script at files/docker-entrypoint.sh)
# COPY files/docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# Prepare dirs and permissions
RUN mkdir -p /docker-entrypoint-initdb.d /var/lib/postgresql/data && chown -R postgres:postgres /var/lib/postgresql /docker-entrypoint-initdb.d /usr/share/postgresql/${PG_MAJOR}/extension

RUN ldconfig

EXPOSE 5432
# adding this in compose file instead
# VOLUME ["/var/lib/postgresql/data"]

# entrypoint script needs root user to init db, then drops to postgres
# USER postgres
WORKDIR /var/lib/postgresql

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["postgres", "-c", "listen_addresses=*", "-D", "/var/lib/postgresql/data"]
