#!/bin/bash

cd eic
curl -L https://github.com/eic/eic-shell/raw/main/install.sh | bash -s -- -c jug_xl -v 23.08.0-stable
