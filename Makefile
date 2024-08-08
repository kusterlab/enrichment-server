IMAGE ?= enrichment_server

INTERACTIVE := $(shell [ -t 0 ] && echo 1)
ifdef INTERACTIVE
	USE_TTY= -t
else
	USE_TTY=
endif

DOCKER_CMD ?= docker run -i ${USE_TTY}

build:
	docker build -f Dockerfile -t $(IMAGE) . || (exit 1)

jump: 
	$(DOCKER_CMD) --entrypoint /bin/bash \
		$(IMAGE)

integration_test:
	$(DOCKER_CMD) --entrypoint "" \
		$(IMAGE) /bin/bash -c "Xvfb :1 -screen 0 1024x768x24 & poetry run python -m pytest -v"
