FROM ubuntu:20.04

# Install system dependencies
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3-dev \
        python3-pip \
        build-essential \
        git \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Copy and install Python dependencies
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy the application code
COPY . .

# Set the command to run the application
CMD [ "python3", "app.py" ]