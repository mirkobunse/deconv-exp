GROUP=`id --group --name`
GID=`id --group`
USER=`id --user --name`
UID=`id --user`

TAG_NAME=$(USER)/mt


.PHONY: default image pngs
default: image


image:
	@ echo "The image '${TAG_NAME}' is built for the Docker image repository '${DOCKER_REPOSITORY}'."
	@ echo "You can set the name of this repository with the environment variable DOCKER_REPOSITORY.\n"
	- docker rmi -f $(TAG_NAME)
	- docker rmi -f $(DOCKER_REPOSITORY)/$(TAG_NAME)
	docker build \
	    --build-arg group=$(GROUP) \
	    --build-arg gid=$(GID) \
	    --build-arg user=$(USER) \
	    --build-arg uid=$(UID) \
	    --tag $(TAG_NAME) \
	    docker/
	docker tag $(TAG_NAME) $(DOCKER_REPOSITORY)/$(TAG_NAME)
	docker push $(DOCKER_REPOSITORY)/$(TAG_NAME)
	docker pull $(DOCKER_REPOSITORY)/$(TAG_NAME)


pngs: $(patsubst res/pdf/%.pdf,res/png/%.png, $(wildcard res/pdf/*.pdf))
res/png/%.png: res/pdf/%.pdf
	convert -density 300 $< $@

