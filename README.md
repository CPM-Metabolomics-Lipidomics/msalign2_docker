# MSalign2

This docker container runs the msalign2 web service.

# Run 

Run the build container with. You need a folder for:

* data
* temporary files

```
docker run -d -p 80:80 -v /your/folder/with/data:/var/www/html/data -v /your/folder/for/temp:/var/www/html/temp ricolumc/msalign2:latest
```

# Usage

Go to: http://localhost/index.php

# Supported files

* mzXML
* mzML (uncompressed or zlib compressed, 32 or 64 bit, indexed or not). MS-Numpress compressed mzML files are not supported.

Master file and files to align can be a mix of mzXML and mzML. The aligned files are written next to the originals:

* mzXML: `<file>.mzXML_ALIGNED.mzxml`
* mzML: `<file>_ALIGNED.mzML` (indexed mzML). The chromatograms (TIC, BPC) are removed, because their time axis is not aligned.
