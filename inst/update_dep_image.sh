# The arm64 image is too large to load with Github actions
# So we build and push locally with the following commands


IMAGE_DEPS=docker.io/breedinginsight/bigapp-deps
DEPS_TAG=r4.5-bioc3.21-2025-08 # Update date here

# amd64 (load locally instead of pushing)
docker buildx build \
  -f Dockerfile.deps \
  --platform linux/amd64 \
  -t $IMAGE_DEPS:$DEPS_TAG-amd64 \
  --progress=plain \
  --load \
  .

# arm64 (push or save to a tar; cannot load multi-arch to local daemon)
docker buildx build \
  -f Dockerfile.deps \
  --platform linux/arm64 \
  -t $IMAGE_DEPS:$DEPS_TAG-arm64 \
  --progress=plain \
  --push \
  .
  
docker buildx imagetools create \
  -t $IMAGE_DEPS:$DEPS_TAG \
  $IMAGE_DEPS:$DEPS_TAG-amd64 \
  $IMAGE_DEPS:$DEPS_TAG-arm64