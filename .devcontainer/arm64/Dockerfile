# Use native multi-arch image (works on ARM64 and AMD64)
# Note: rocker/shiny-verse has limited ARM64 support because Shiny Server
# doesn't distribute ARM64 binaries. Use rocker/r-ver instead.
FROM rocker/r-ver:4.4.3

# 1. Configure Posit Package Manager for prebuilt ARM64 binaries
# This dramatically speeds up package installation (minutes instead of 30+ min)
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))' >> /usr/local/lib/R/etc/Rprofile.site

# 2. System deps (still needed for some packages that compile from source)
RUN apt-get update && apt-get install -y \
    build-essential \
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
    cmake \
    git \
    curl \
 && rm -rf /var/lib/apt/lists/*

# 3. Install Docker CLI (for spawning sibling MCP containers)
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg && \
    echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu jammy stable" > /etc/apt/sources.list.d/docker.list && \
    apt-get update && apt-get install -y docker-ce-cli && \
    rm -rf /var/lib/apt/lists/*

# 4. Install all R packages needed for analysis, IDE integration, and MCP
# Note: dependencies=NA skips optional "Suggests" (e.g., otelsdk) for faster builds
RUN R -e "options(warn = 2); install.packages(c('tidyverse', 'shiny', 'readxl', 'fst', 'jsonlite', 'plotly', 'randomForest', 'nnet', 'caret', 'pROC', 'dbscan', 'umap', 'igraph', 'fastDummies', 'smotefamily', 'rstudioapi', 'languageserver', 'btw', 'ellmer', 'mcptools', 'Rprofet'), dependencies = NA)"

# 7. Install Node.js LTS from NodeSource
RUN curl -fsSL https://deb.nodesource.com/setup_lts.x | bash - \
    && apt-get install -y nodejs \
    && npm install -g npm@latest

# 8. Install Claude Code globally
RUN npm install -g @anthropic-ai/claude-code
RUN npm install -g @openai/codex

# 9. Additional R packages can be installed through a packages.R file
#COPY ./packages.R /tmp/packages.R
#RUN Rscript /tmp/packages.R

# 10. Give shiny passwordless sudo for updating Claude Code from VS Code
RUN useradd -m -s /bin/bash shiny && \
    echo 'shiny ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers && \
    usermod -aG root shiny && \
    usermod -aG staff shiny && \
    rm -rf /var/lib/apt/lists/*

# Expose Shiny server port
EXPOSE 3838
