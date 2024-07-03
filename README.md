# GenomeProt: an integrated proteogenomics data analysis platform for long-read RNA-Seq datasets

### Obtaining the shiny application with Docker
Make sure you have Docker (link) installed and the application running in the background before you beign.

Open your terminal application and run:
```
docker image pull josiegleeson/GenomeProt:latest
```

### Running the shiny application with Docker

Open your terminal application and run:
```
docker run --rm -p 3838:3838 GenomeProt
```
The --rm removes the container after it’s stopped and the -p 3838:3838 maps your local port 3838, to the same port inside the container.

To **access the local shiny application**, navigate to this link on your web browser http://0.0.0.0:3838.

You can now upload all files and run the steps in your web browser.

To stop the container, head back to the terminal where docker is running and press ctrl+c.
