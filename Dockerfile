# syntax=docker/dockerfile:1
# check=skip=SecretsUsedInArgOrEnv,JSONArgsRecommended

FROM anaconda/miniconda:latest

# is_install_fragpipe: to install FragPipe, set this to "true"; any other value means FragPipe will not be installed
ARG is_install_fragpipe=true

# fragpipe_token: the 6-digit token used to install FragPipe tools
ARG fragpipe_token=123456

ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true

# Install gcc because it's necessary for building biopython
# Install unzip and gzip for unpacking the reference and test datasets

RUN apt update && \
    apt upgrade && \
    apt install --yes gcc unzip gzip && \
    useradd --create-home --shell /bin/bash user

# Set up the user directory

RUN mkdir -p /home/user/Desktop/GenomeProt
WORKDIR /home/user/Desktop/GenomeProt
COPY . .
RUN chown -R user:user GenomeProt

# Use the new user account from here onwards

USER user

# Install FragPipe if necessary

RUN if [ "$is_install_fragpipe" = "true" ]; then \
        printf "Y\nY\n$fragpipe_token\n" | python3 fragpipe_installer.py /home/user/Desktop/GenomeProt; \
    fi

# Build the conda environment

RUN if [ "$is_install_fragpipe" = "true" ]; then \
        conda env create -f conda_env.yaml; \
    else \
        conda env create -f conda_env_no_fragpipe.yaml; \
    fi

# Initialize the conda environment and activate it in the user's shell

RUN conda init && \
    echo 'conda activate GenomeProt_env' >> ~/.bashrc

# Prepare the reference and test datasets

RUN bash prepare_ref_and_test_data.sh

# Expose the port GenomeProt will run on to the host system

EXPOSE 3838

# Run GenomeProt

CMD /bin/bash -c 'conda run -n GenomeProt_env Rscript -e "shiny::runApp('"'"'GenomeProt/'"'"', host='"'"'0.0.0.0'"'"', port=3838)"'