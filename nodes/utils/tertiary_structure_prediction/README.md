# Setup tertiary structure prediction workflow

In order to run the Psi-Blast-based tertiary structure prediction workflow (TSPW), a few things 
needs to be set up: RaptorX, Modeller and SpineX. Modeller is needed for the comparative modelling of
RaptorX. SpineX will be used to generate secondary structure information from the predicted tertiary structure.
Moreover, RaptorX generates predictions from known structure, such that some databases are also required. This
document guides you through all necessary steps.

## RaptorX

Setup RaptorX, including the databases. A comprehensive manual can be found in the download section of
RaptorX (see section `Registration`).

### Registration

First, register under http://raptorx.uchicago.edu/download/ to obtain your personal download link.

### Databases

The following databases are required. Download, unzip and save these in the respective location.

| Name | Link  | Directory  |
|---|---|---|
|  a | a  | a  |
|   |   |   |
|   |   |   |

### Initialize

1) Download the actual tool (CNFsearch), unzip it and save it as `peptidereactor/RaptorX/`.
2) `cd` into `peptidereactor/RaptorX/` and run `setup.pl`. 

### RaptorX

1) Register under http://raptorx.uchicago.edu/download/. 
2) Paste link into file, e.g., `download_links/raptorx_download_link.txt`. 
3) Access it via `config["raptorx_download_link_in"]`.

### Modeller

1) Obtain the license key at https://salilab.org/modeller/registration.html.
2) Paste key into a file, e.g., `download_links/modeller_license_key.txt`.
3) Access it via `config["modeller_license_key_in"]`.