VERSION = 1.2
TAG1 = quay.io/ymostovoy/lr-palmer-filter:$(VERSION)
TAG2 = quay.io/ymostovoy/lr-palmer-filter:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
