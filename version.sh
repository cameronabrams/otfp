#!/bin/sh

# Get version
version=$(git describe --abbrev=4 --dirty --always 2>/dev/null || cat version 2>/dev/null || echo unknown)

# Update version file
[ -z "${MESON_DIST_ROOT:-}" ] || echo "$version" > "${MESON_DIST_ROOT}/version"

# Print version
echo $version
