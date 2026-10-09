[//]: # ([![DOI]&#40;https://zenodo.org/badge/DOI/10.5281/zenodo.17650102.svg&#41;]&#40;https://doi.org/10.5281/zenodo.17650102&#41;)

![GitHub License](https://img.shields.io/github/license/sophietheis/KineticAnalysis)

[//]: # ([![PyPI]&#40;https://img.shields.io/pypi/v/epitools.svg?color=green&#41;]&#40;https://pypi.org/project/epitools&#41;)

[//]: # ([![tests]&#40;https://github.com/epitools/epitools/actions/workflows/test.yml/badge.svg&#41;]&#40;https://github.com/epitools/epitools/actions/workflows/test.yml&#41;)

[//]: # ([![Documentation]&#40;https://readthedocs.org/projects/epitools/badge/?version=latest&#41;]&#40;https://epitools.readthedocs.io/en/latest/?badge=latest&#41;)

_______
[![Docker](https://img.shields.io/badge/Docker-2496ED?logo=docker&logoColor=fff)](https://hub.docker.com/repository/docker/sophiets/kineticapp/general)
![Docker Image Version](https://img.shields.io/docker/v/sophiets/kineticapp)
![Docker Pulls](https://img.shields.io/docker/pulls/sophiets/kineticapp)


_____

# Overview

`Kinetic analysis` is an interactive web application for analysing mRNA translation dynamics using the **SunTag fluorescence system**. Built with [Plotly Dash](https://dash.plotly.com/), it provides a complete, researcher-facing analysis pipeline — from experimental design and synthetic data generation, through kinetic parameter fitting, ribosome counting, and diffusion analysis, to population-level result aggregation.

`Kinetic analysis` is designed for research scientists studying translation kinetics via live-cell single-molecule fluorescence imaging. No programming experience is required — all analysis is performed through a tab-based graphical interface.


This project is a collaboration between [Mounia Lagha's](http://www.laghalab.com/) and [Tim Saunders'](https://mechanochemistry.org/Saunders/MainSite/Saunders_lab_v4_3.htm) teams. 


# Features
| Module | Description |
|---|---|
| **Acquisition Parameters** | Calculate optimal imaging time step, recording duration for your experimental setup |
| **Generate Tracks** | Simulate synthetic polysome fluorescence tracks for pipeline benchmarking |
| **Equation Choice** | Visualise and compare autocorrelation equation term contributions to select the best fitting model |
| **Track Analysis** | Fit kinetic parameters (elongation rate *k*_elong, initiation rate *c*_init) to experimental fluorescence tracks via autocorrelation |
| **Combine Results** | Merge and aggregate tracks |
| **Count Ribosomes** | Normalise fluorescence intensity to estimate ribosome counts per frame |
| **MSD Analysis** | Compute Mean Squared Displacement and extract diffusion coefficients from trajectory data |



# Getting started

## 🐳 Docker

The easiest way to run `kinetic analysis` is via the pre-built Docker image published on Docker Hub — no Python installation, no dependency management required.

<details>
<summary>Open to see how to install </summary>

### Prerequisites

Install **Docker Desktop** for your operating system:
- [Download for Windows](https://docs.docker.com/desktop/install/windows-install/)
- [Download for macOS](https://docs.docker.com/desktop/install/mac-install/)
- [Download for Linux](https://docs.docker.com/desktop/install/linux-install/)

Once installed, open Docker Desktop and make sure it is running (the whale icon appears in your taskbar/menu bar).

### Step 1 — Pull the `kineticapp` image

1. Open **Docker Desktop**
2. Click the **Search** bar at the top
3. Type `sophiets/kineticapp` and press Enter
4. Find the image `sophiets/kineticapp` in the results
5. Click **Pull**

Docker will download the image automatically. You will see a progress bar. Wait until it says **Pull complete**.

> You can also find the image on Docker Hub:
> 👉 [hub.docker.com/r/sophiets/kineticapp](https://hub.docker.com/r/sophiets/kineticapp)

### Step 2 — Run the container

1. Go to the **Images** tab in the left sidebar
2. Find `sophiets/kineticapp` in the list
3. Click the ▶️ **Run** button on the right

A configuration panel will appear — click **Optional settings** to expand it:

| Setting | Value to enter |
|---|---|
| **Container name** | `kinetic-app` *(or any name you like)* |
| **Host port** | `5001` |

Leave everything else as default. Click **Run**.

### Step 3 — Open the app

Once the container is running, open your web browser and go to:

```
http://localhost:5001
```

`KineticApp` will load in your browser.

### Step 4 — Stop the app

1. Go to the **Containers** tab in the left sidebar
2. Find `photon-app` in the list
3. Click the ⏹️ **Stop** button

To start it again later, click the ▶️ **Start** button — no need to pull or configure again.


### Step 5 — Update to the latest version

When a new version of PHOTON is released:

1. Go to the **Images** tab
2. Find `sophiets/kineticapp`
3. Click the ⋮ **menu** → **Pull** to download the latest version
4. Delete the old container in the **Containers** tab and run a new one from the updated image (repeat Step 2)

### Troubleshooting

| Problem | Solution |
|---|---|
| Port 5001 already in use | Change the host port to `5002` or any free port, then open `http://localhost:5002` |
| App not loading in browser | Make sure the container status shows **Running** (green) in Docker Desktop |
| Docker Desktop not starting | Ensure virtualisation is enabled in your BIOS / system settings |

</details>
<!-- Docker is the recommended way to run `kinetic analysis` in a reproducible, dependency-free environment — no Python installation required. -->


<!-- See [docker install](https://github.com/sophietheis/KineticAnalysis/blob/main/doc/Install_docker_kineticapp.pdf) -->

## Python method

<details>
<summary>Open to see how to install </summary>
Clone the repository.
```sh
git clone https://github.com/sophietheis/KineticAnalysis
```

We recommend to create an environment to install `KineticAnalysis` package. 

```sh
conda create --name kinetic-env python=3.11
conda activate kinetic-env
```

<!-- The recommended way to install `KineticAnalysis` is via [pip](https://pypi.org/project/pip/). -->


To install all dependencies:
```sh
pip install --no-cache-dir -r requirements.txt
```

To install the latest development version, clone this repository and run:
```sh
pip install .
```

Running the app
```sh
python3 src/kinetic_analysis/app_dash.py 
```
Open your favorite web browser and enter the address: http://127.0.0.1:5001/

</details>

## Issues

If you encounter any problems, please [open an issue](https://github.com/sophietheis/KineticAnalysis/issues) along with detailed descriptions.

## Contributing
Contributions are very welcome.


