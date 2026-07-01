FROM quay.io/jupyter/minimal-notebook:2025-04-14
LABEL authors="Max Winokan"


# Ported from Max's Dockerfile to support local development

# combines conda and pip environments. There's a conflict somewhere
# which I haven't resolved (TODO) between something in conda and pip
# envs, so cannot use a single env

# this is what I'm stuck on. Project depends on numpy v1.23 (why?) and
# I can't go over 3.12
ARG PYTHON_VERSION=3.12


WORKDIR "/home/code/HIPPO"
# WORKDIR "/home/code/HIPPO/data"
# WORKDIR "/home/code/HIPPO/hippo"

USER 0


# install package dependencies into different virtual env
COPY uv.lock pyproject.toml ./
RUN pip install --upgrade pip && python -m pip install uv

# install all dependencies into active environment without updating lockfile
# RUN uv sync --frozen --active
RUN uv venv /home/code/HIPPO/.venv
RUN uv sync --frozen --python /home/code/HIPPO/.venv/bin/python

# now add venv python to path so conda python can find it
ENV PATH="/home/code/HIPPO/.venv/bin:$PATH"
ENV PYTHONPATH="/home/code/HIPPO/.venv/lib/python${PYTHON_VERSION}/site-packages:$PYTHONPATH"


# patch rich
RUN python -c "import mrich; mrich.patch_rich_jupyter_margins()"


# notebooks
RUN chown ${NB_USER} "/home/code" && sudo apt update && sudo apt install screen -y

# force numpy update. since there's now a new venv, make sure conda is
# used (othewise numpy will be invisible)
RUN uv pip install --upgrade numpy==2.4.4
USER ${NB_USER}


# WORKDIR "/home/code/HIPPO"
