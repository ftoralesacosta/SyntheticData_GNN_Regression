#!/bin/bash
cd eic
curl https://eicweb.phy.anl.gov/containers/eic_container/-/raw/master/install.sh | bash
mkdir development

./eic-shell
