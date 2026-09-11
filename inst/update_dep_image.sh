#!/bin/bash
set -euo pipefail

# The arm64 image is too large to load with Github actions
# So we build and push locally with the following commands

# There is a file inside the image with all packages versions listed
# find it at /opt/pins/deps_versions.tsv

# Make sure to login with BI user
IMAGE_DEPS=docker.io/breedinginsight/bigapp-deps
BUILD_ID=build-$(date +%s) # throwaway tag used only to inspect the freshly built image

# amd64 (load locally so we can inspect it before deciding the real tag)
docker buildx build \
  -f Dockerfile.deps \
  --platform linux/amd64 \
  -t $IMAGE_DEPS:$BUILD_ID-amd64 \
  --progress=plain \
  --load \
  .

# Derive DEPS_TAG from what's actually inside the image (R version, Bioc version, build date)
# instead of hardcoding it, so the tag can never drift from reality.
R_VER=$(docker run --rm $IMAGE_DEPS:$BUILD_ID-amd64 Rscript -e 'cat(R.version$major, sub("\\..*", "", R.version$minor), sep=".")')
BIOC_VER=$(docker run --rm $IMAGE_DEPS:$BUILD_ID-amd64 Rscript -e 'cat(as.character(BiocManager::version()))')
DEPS_TAG="r${R_VER}-bioc${BIOC_VER}-$(date +%Y-%m)"
echo "Resolved DEPS_TAG=$DEPS_TAG"

docker tag $IMAGE_DEPS:$BUILD_ID-amd64 $IMAGE_DEPS:$DEPS_TAG-amd64
docker rmi $IMAGE_DEPS:$BUILD_ID-amd64
docker push $IMAGE_DEPS:$DEPS_TAG-amd64

# arm64 (push or save to a tar; cannot load multi-arch to local daemon)
docker buildx build \
  -f Dockerfile.deps \
  --platform linux/arm64 \
  --build-arg DEPS_TAG=$DEPS_TAG \
  -t $IMAGE_DEPS:$DEPS_TAG-arm64 \
  --progress=plain \
  --push \
  .

docker buildx imagetools create \
  -t $IMAGE_DEPS:$DEPS_TAG \
  -t $IMAGE_DEPS:latest \
  $IMAGE_DEPS:$DEPS_TAG-amd64 \
  $IMAGE_DEPS:$DEPS_TAG-arm64

echo "Pushed $IMAGE_DEPS:$DEPS_TAG and $IMAGE_DEPS:latest"