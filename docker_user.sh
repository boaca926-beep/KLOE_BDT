#!/bin/bash
#Add user to docker group, avoid sudo
#sudo usermod -aG docker $USER

# Log out and back in, or run:
#newgrp docker

# Verify it works
docker ps
