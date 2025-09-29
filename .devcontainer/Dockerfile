# Choose your R version (update as you like)
FROM gexijin/idep:latest

# System deps commonly needed by R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev libssl-dev libxml2-dev git curl && \
    rm -rf /var/lib/apt/lists/*

# R packages for VS Code integration: language server + debugger
RUN R -q -e 'install.packages("vscDebugger", repos=c("https://manuelhentschel.r-universe.dev","https://cloud.r-project.org"))'
RUN R -q -e 'install.packages("languageserver", repos="https://cloud.r-project.org")'

# 2. Install Node.js LTS from NodeSource
RUN curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g npm@latest
    
# Install Claude Code globally
RUN npm install -g @anthropic-ai/claude-code

# Install additional R packages for the script
#COPY ./install_packages.R /tmp/install_packages.R
#RUN Rscript /tmp/install_packages.R

EXPOSE 3838