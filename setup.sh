#!/bin/bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
curl -o clustalo 'https://web.archive.org/web/20120110042415/http://www.clustal.org/omega/clustalo-1.0.3-Ubuntu-amd64'
chmod +x datasets dataformat clustalo
echo "Setup complete."

