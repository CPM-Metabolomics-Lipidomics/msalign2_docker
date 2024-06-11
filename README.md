# MSalign2

This docker container runs the msalign2 web service.

# Run 

Run the build container with. You need a folder for:

* data
* temporary files

```
docker run -d -p 80:80 -v /your/folder/with/data:/var/www/html/data -v /your/folder/for/temp:/var/www/html/temp ricolumc/msalign2:latest
```
