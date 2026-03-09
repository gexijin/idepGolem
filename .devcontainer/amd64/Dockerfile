# choose a Dockerhub base image with R, Shiny Server, and tidyverse packages
FROM rocker/shiny-verse:latest

# 1. System deps commonly needed by R packages (including cmake for nanonext/mcptools)
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libsodium-dev \
    libcairo2-dev \
    git \
    curl \
    gnupg \
 && rm -rf /var/lib/apt/lists/*

# 2. Install Docker CLI (for spawning sibling MCP containers)
# Note: rocker/shiny-verse is Ubuntu-based, so we use the Ubuntu Docker repo
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu $(. /etc/os-release && echo "$VERSION_CODENAME") stable" > /etc/apt/sources.list.d/docker.list && \
    apt-get update && apt-get install -y docker-ce-cli && \
    rm -rf /var/lib/apt/lists/*

## 3. Install all R packages needed for analysis, IDE integration, and MCP
RUN R -q -e 'install.packages(c("readxl", "fst", "jsonlite", "plotly", "randomForest", "nnet", "caret", "pROC", "dbscan", "umap", "igraph", "fastDummies", "smotefamily", "rstudioapi", "languageserver", "btw", "ellmer", "mcptools", "Rprofet"), repos="https://cloud.r-project.org")'

# 5. Install Node.js LTS from NodeSource
RUN curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g npm@latest

# 6. Install Claude Code globally
RUN npm install -g @anthropic-ai/claude-code
RUN npm install -g @openai/codex

# 7. Additional R packages can be installed through a packages.R file
#COPY ./packages.R /tmp/packages.R
#RUN Rscript /tmp/packages.R

# 8. Give shiny passwordless sudo and Docker socket access
# Add shiny to root group for Docker socket access (socket is typically root:root or root:docker)
RUN echo 'shiny ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers && \
    usermod -aG root shiny && \
    usermod -aG staff shiny && \
    rm -rf /var/lib/apt/lists/*

# Expose Shiny server port
EXPOSE 3838
