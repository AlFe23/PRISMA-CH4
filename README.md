 

# PRISMA-CH4

This repository contains the code and resources for detecting and quantifying methane (CH<sub>4</sub>) emissions using PRISMA hyperspectral satellite imagery.

The methodology, originally developed for the CLEAR-UP project funded by ASI, leverages the PRISMA Hyperspectral sensor to monitor methane emissions using an enhanced matched filter technique known as the Cluster Tuned Matched Filter (CTMF). This is applied to L1 radiance images. The unit absorption spectrum, obtained using multiple MODTRAN6® simulations, serves as the target in the matched filter, enabling methane detection and quantification. Key steps include utilizing PRISMA L1 radiance images in the 2300 nm short-wave infrared methane absorption window, adapting the matched filter to account for background clutter, and applying clustering for refined background spectrum calculations. Methane concentration is estimated by linearizing the Beer-Lambert absorption law and using the unit absorption spectrum, which represents the unit methane absorption as a function of wavelength. Automation through a lookup table (LUT) based on precomputed radiances for various conditions enables real-time applicability. Results computed over 200 PRISMA acquisitions revealed methane plumes at various sites of interest in countries including Algeria, Argentina, Brazil, China, India, Mexico, Turkmenistan, and the USA. For instance, the Natural Gas compression plant in Kamyshlydzha, Turkmenistan, consistently exhibited methane plumes across all 13 images analyzed. Several large landfills were examined using available PRISMA images or through specific acquisition requests. Evidence of methane concentration enhancement was observed in the Buenos Aires landfill and the Pirana landfill in Ahmedabad.

Check also the implementation for EnMAP: [EnMAP-CH4](https://github.com/AlFe23/EnMAP-CH4)

<div align="center">
<img src="https://github.com/AlFe23/PRISMA-CH4/assets/105355911/87a5817e-41c2-4d25-9496-8b8e5d2e71fe" width="20%">
</div>

**Authors:** Alvise Ferrari, Giovanni Laneve 

**Contact:** alvise.ferrari@uniroma1.it, ferrarialvise@gmail.com, giovanni.laneve@uniroma1.it

**Citation:**

- A. Ferrari, G. Laneve, V. Pampanoni, A. Carvajal and F. Rossi, "Monitoring Methane Emissions from Landfills Using Prisma Imagery," IGARSS 2024 - 2024 IEEE International Geoscience and Remote Sensing Symposium, Athens, Greece, 2024, pp. 3663-3667, doi: 10.1109/IGARSS53475.2024.10642079.

- Ferrari, A., Laneve, G., Tellez, R.A.C., Pampanoni, V., Saquella, S. and Guarini, R., 2025. Analysis of Local Methane Emissions Using Near-Simultaneous Multi-Satellite Observations: Insights from Landfills and Oil-Gas Facilities. arXiv preprint arXiv:2506.01113.

**Links:**

- [Docker prisma-ch4-mapper.tar](https://drive.google.com/file/d/1xFfn2pqW5srqmd0dexAtlwD8hRF18fMS/view?usp=sharing)
- [LUT CH4](https://drive.google.com/file/d/196adGp_XCcTXAk3SRjiOnBJxUhDANNvn/view?usp=sharing)
- [DEM](https://drive.google.com/file/d/10e1VtibryVxcHT4-Gb0ryhyk17JCF04f/view?usp=sharing)
- [PRISMA Sample Images](https://drive.google.com/drive/folders/1IqJE_szLeWtHDRRdjURAl2JIBeQL_9zy?usp=sharing)


## Index

1. **Algorithm Theoretical Basis**  
   - 1.1 Clutter Matched Filter  
   - 1.2 Concentration Estimation  
   - 1.3 Scene-Specific Target Spectrum Automatic Generation  
   - 1.4 Processing Flow Chart  
   - 1.5 Examples of Automatically Generated Outputs  

2. **User Manual for prisma-ch4-mapper Docker Container (v2)**  
   - 2.1 Provided Files  
   - 2.2 Docker Container Usage  
   - 2.3 Example Docker Run Commands  
     - 2.3.1 Using `run_ctmf_prisma.py`  
     - 2.3.2 Using plain `docker`  
   - 2.4 Input Definition  
   - 2.5 Output Definition  
   - 2.6 Future Updates  

3. **References**


## 1. **Algorithm Theoretical Basis**

### 1.1 Clutter Matched Filter

In the SWIR spectrum, the absorption windows of methane are centered about 1700 nm and 2300 nm. The radiance measured by a hyperspectral sensor is modeled as a superposition of signal $\( b \)$, scaled by its strength $\( \alpha \)$, an average background radiance $\( L_{\text{mean}} \)$ (averaged over the entire image), and a zero-mean noise or clutter term $\( \epsilon \)$ [12][19].

<div align="center">
  
$`r = \alpha b + L_{\text{mean}} + \epsilon`$

</div>

$\( \epsilon \)$ represents both sensor noise and scene clutter (non-desirable components). For an uncorrelated background, the optimal filter is the basic matched filter (i.e., the target itself, scaled to make the variance of the filtered image equal to one). Considering the real case of background clutter with correlation between spectral channels, the optimal filter is “matched” both to the signature of the target and the background:

**Clutter Matched Filter (CMF):**

<div align="center">
 
$`q = \frac{C^{-1} b}{\sqrt{b^T C^{-1} b}}`$

</div>

Here, $\( q \)$ is normalized to ensure that in the absence of signal, the variance of the matched filter image, $\( q^T r_i \)$ (for the i-th pixel) is one. Values larger than one indicate strong evidence of the signature's presence, quantified as a "number of sigmas", assuming that the clutter itself is Gaussian. When the clutter's covariance matrix $\( C \)$ is accurately known, the equation defines the optimal matched filter, maximizing the signal-to-clutter ratio. However, this covariance is seldom known beforehand and is typically deduced from the data, often by averaging the outer product of the mean-subtracted radiance across all pixels.

<div align="center">

$`C \approx \frac{1}{N} \sum_{i=1}^N (L_i - L_{\text{mean}})(L_i - L_{\text{mean}})^T`$

</div>

**Cluster Tuned Matched Filter (CTMF):**

The CTMF [12], performs image clustering before applying the matched filter [14][15][16][17][18]. The average background spectrum and covariance matrix are calculated for each class. The CTMF score for the pixel $\((x,y)\)$ in the j-th class is given by:

<div align="center">

$`\alpha_j (x,y) = q_j (L_{(x,y)_j} - L_{\text{mean}_j})`$

</div>

<div align="center">

$`q_j = \frac{b_j^T C_j^{-1}}{\sqrt{b_j^T C_j^{-1} b_j}}`$

</div>

Where the optimal filter $\( q_j \)$ for the j-th class is expressed as [12], $\( L_{\text{mean},j} \)$ is the average radiance of class $j$ and $\( C_j \)$ is the related covariance matrix [18] and $\( b_j \)$ is the target signal. An important assumption is made approximating the off-plume radiance with the average radiance of class $j$, $\( L_{\text{mean},j} \)$, to determine the expression of $\( b \)$.

<div align="center">

$`L_{\text{no\_pl}} = L_{\text{mean},j}`$

</div>

### 1.2 Concentration Estimation (v2)

The methane column enhancement $`\Delta X_{\mathrm{CH}_4}`$ (in ppm·m) is obtained by linearizing the Beer–Lambert law for optically thin plumes ($`\tau \lesssim 0.05`$) and applying a maximum-likelihood matched filter:

<div align="center">

$`\Delta X_{\mathrm{CH}_4}(x,y) =
\frac{(L_{(x,y)} - \mu)^\top\,\Sigma^{-1}\,t}
     {t^\top\,\Sigma^{-1}\,t}`$

</div>

<div align="center">

$`t = -A \bigl(1 + \sec\theta\bigr)\,\mu`$

</div>

where  
* **$`L_{(x,y)}`$** is the at-sensor radiance vector,  
* **$`\mu, \Sigma`$** are the background mean and covariance (global or per-cluster),  
* **$`t`$** is the unit-target spectrum ($`\mathrm{radiance}\cdot(\mathrm{ppm}\,\mathrm{m})^{-1}`$) from LUT/MODTRAN,  
* **$`A`$** is the scene-specific $`\mathrm{CH}_4`$ absorption coefficient,  
* **$`\theta`$** is the solar zenith angle.  

Because $`\Sigma^{-1}`$ appears in both numerator and denominator, the result is already in physical units (ppm·m) without empirical scaling.  

Use `--n_clusters k` to tune $`\mu`$, $`\Sigma`$ and $`t`$ per cluster for improved SNR on heterogeneous scenes.  



### 1.3 Scene-Specific Target Spectrum Automatic Generation

The unit absorption spectrum is derived from multiple simulations with varying methane gas concentrations using radiative transfer models such as MODTRAN6®. This spectrum represents the methane absorption as a function of wavelength, normalized per unit path length and concentration. The unit for this measurement is typically $(ppm·m)\(^{-1}\)$, indicating the inverse of the product of methane concentration (in parts per million) and the path length (in meters) over which the absorption occurs.

Subsequently, the simulated radiance spectra are convolved with the PRISMA spectral response function, and a regression analysis is carried out to obtain the unit absorption spectrum for each PRISMA spectral channel. The slope of the regression line, which best represents the relationship between the concentration-path product and the natural logarithm of the radiance, provides the unit absorption value for each specific wavelength band. This unit absorption spectrum is scaled by element-wise multiplication with the average radiance at each wavelength to generate the target spectrum, which is then used in the matched filter. To automate the generation of the unit target spectrum and make it independent of real-time MODTRAN6® runs, a lookup table (LUT) has been developed based on the methodology described by Foote et al. (2020) [4]. This LUT contains precomputed at-sensor radiances for a wide range of atmospheric and geometric parameters, including variations in sensor altitude, water vapor content, ground elevation, sun zenith angle, and methane concentration. For each PRISMA acquisition, the Level 2C PRISMA image is also used to read the water vapor content, and a Digital Elevation Model (DEM) is used to read the ground elevation. By interpolating within the LUT, it is possible to quickly generate a unit target spectrum that accurately matches the conditions of a given satellite pass.

<div align="center">
<img src="https://github.com/AlFe23/PRISMA-CH4/assets/105355911/160778eb-f03e-477b-83be-0ce0200c0bbb" width="50%">
</div>


### 1.4 Processing Flow Chart

Below a figure reporting the processing workflow implemented:

<div align="center">
<img src="https://github.com/AlFe23/PRISMA-CH4/assets/105355911/631d12fe-f2eb-424c-9e27-da3ace35469a" width="60%">
</div>

### 1.5 Examples of Automatically Generated Outputs

<div align="center">
<img src="https://github.com/AlFe23/PRISMA-CH4/assets/105355911/d471e0b0-3343-4c26-a5ef-7cb7f16bfa4d" width="50%">
</div>




##  2. **User Manual for `prisma-ch4-mapper` Docker Container**

### 2.1 Provided Files:

- `prisma-ch4-mapper.tar` (Docker image)
- `dataset_ch4_full.hdf5`
- `srtm30plus_v11_land.nc`
- sample PRISMA L1 and L2C images in areas with potential CH4 emitters
- Python script `run_ctmf_prisma.py` for ease the docker container use


All files can be found at this link, navigating to the folder `2.Rilevamento CH4`.

### 2.2 Docker Container Usage

The Docker image (`prisma-ch4-mapper.tar`) now supports two modes of operation:

- **Single mode**: process one acquisition at a time by supplying the paths to the extracted `.he5` L1 and L2C files.  
- **Batch mode**: process a directory tree of PRISMA acquisitions (each in its own folder, either zipped or already extracted) in one go.

You can drive the container either via the helper script `run_ctmf_prisma.py` (recommended) or with plain `docker` commands.

---

### 2.3 Example Docker Run Commands

#### 2.3.1 Using `run_ctmf_prisma.py`

> **Note:** place `prisma-ch4-mapper.tar` in the same folder as `run_ctmf_prisma.py`.

##### Single mode

```bash
python3 run_ctmf_prisma.py \
  -t prisma-ch4-mapper.tar \
  single \
    /path/to/PRS_L1_STD_OFFL_xxx.he5 \
    /path/to/PRS_L2C_STD_xxx.he5 \
    /path/to/srtm30plus_v11_land.nc \
    /path/to/dataset_ch4_full.hdf5 \
    /path/to/output/folder \
    -k 1 --min_wl 2100 --max_wl 2450
```

##### Batch mode

```bash
python3 run_ctmf_prisma.py \
  -t prisma-ch4-mapper.tar \
  batch \
    /path/to/parent_directory_of_acquisitions \
    /path/to/srtm30plus_v11_land.nc \
    /path/to/dataset_ch4_full.hdf5 \
    --output_root /path/to/output/folder \
    -k 1 --min_wl 2100 --max_wl 2450
```

---

#### 2.3.2 Using plain Docker

First load the image:

```bash
docker load -i prisma-ch4-mapper.tar
```

##### Single mode

```bash
docker run --rm \
  -v /host/path/to/data:/data \
  prisma-ch4-mapper:latest single \
    "/data/PRS_L1_STD_OFFL_xxx.he5" \
    "/data/PRS_L2C_STD_xxx.he5" \
    "/data/srtm30plus_v11_land.nc" \
    "/data/dataset_ch4_full.hdf5" \
    "/data/output_folder" \
    -k 1 --min_wl 2100 --max_wl 2450
```

##### Batch mode

```bash
docker run --rm \
  -v /host/path/to/acquisitions:/data/in \
  -v /host/path/to/resources:/data/res \
  -v /host/path/to/output:/data/out \
  prisma-ch4-mapper:latest batch \
    "/data/in" \
    "/data/res/DEM_1Km/srtm30plus_v11_land.nc" \
    "/data/res/LUTs/dataset_ch4_full.hdf5" \
    --output_root "/data/out" \
    -k 1 --min_wl 2100 --max_wl 2450
```

---

> **Spectral window options**  
> You can control the wavelength range with `--min_wl` and `--max_wl` (in nm).  
> **clustering option**  
> You can choose to activate the CTMF (cluster tuned matched filter) simply setting the clustering flag `-k` >1.  
>  
> **Single mode** requires explicit paths to the extracted `.he5` files.  
> **Batch mode** will scan subfolders of the input directory for PRISMA L1 and L2C files (zipped or unzipped) and process them in sequence.  
>  


### 2.3.4 Input Definition

- **L1_file:** Path to the PRISMA L1 file in original format *.he5, as extracted from the zip provided by ASI.
- **L2C_file:** Path to the PRISMA L2C file in original format *.he5, as extracted from the zip provided by ASI.
- **dem_file:** Path to the `srtm30plus_v11_land.nc` file, as provided by EOSIAL.
- **lut_file:** Path to the `dataset_ch4_full.hdf5` file, as provided by EOSIAL.
- **output_dir:** Path to the directory where you want to save the outputs.

### 2.3.5 Output Definition

The output files generated by the software are as follows:

- `PRS_L1_STD_OFFL_XXXXXXXXXXXXXX_XXXXXXXXXXXXXX_0001_rgb.tif`: RGB image of the PRISMA L1 scene
- `PRS_L1_STD_OFFL_XXXXXXXXXXXXXX_XXXXXXXXXXXXXX_0001_classified.tif`: Unsupervised classification image of the scene (k-means)
- `PRS_L1_STD_OFFL_XXXXXXXXXXXXXX_XXXXXXXXXXXXXX_0001_MF.tif`: Matched filter score image, derived from the PRISMA L1 image, indicating the probability of methane presence.
- `PRS_L1_STD_OFFL_XXXXXXXXXXXXXX_XXXXXXXXXXXXXX_0001_MF_concentration.tif`: CH4 concentration map in physical units [ppm∙m]. This corresponds to the CLEAR-UP product 'PRISMA_L1_CH4enhancement' as defined in the official project documentation.

### 2.3.6 Future Updates

Future versions of the script will not require the use of additional input files for generating the 'PRISMA_L1_CH4enhancement' product. It is planned to use a Sentinel-2 image to coregister the final 'PRISMA_L2_CH4enhancement' product, which will be automatically downloaded unless explicitly provided as an additional input. All future released versions will have the same usage structure. Simply download the updated Docker image, without modifying the integration of the Docker within the automatic product generation system or through UI and web apps.


### 3. References

<small>
[1] Saunois et al., "The global methane budget 2000–2012," Earth Syst. Sci. Data, vol. 8, no. 2, pp. 697–751, Dec. 2016, doi: 10.5194/essd-8-697-2016.  
  
[2] Myhre et al., "Anthropogenic and natural radiative forcing," in Climate Change 2013—The Physical Science Basis Working Group I Contribution to the Fifth Assessment Report of the Intergovernmental Panel on Climate Change. Cambridge: Cambridge, U.K.: Cambridge Univ. Press, 2014, pp. 659–740, doi: 10.1017/CBO9781107415324.018.  

[3] Foote, Markus D., et al. "Fast and accurate retrieval of methane concentration from imaging spectrometer data using sparsity prior." IEEE Transactions on Geoscience and Remote Sensing 58.9 (2020): 6480-6492.  

[4] Foote, Markus D., et al. "Impact of scene-specific enhancement spectra on matched filter greenhouse gas retrievals from imaging spectroscopy." Remote Sensing of Environment 264 (2021): 112574.  

[5] Ayasse, Alana K., et al. "Methane mapping with future satellite imaging spectrometers." Remote Sensing 11.24 (2019): 3054. 

[6] Krautwurst, Sven, et al. "Methane emissions from a Californian landfill, determined from airborne remote sensing and in situ measurements." Atmospheric Measurement Techniques 10.9 (2017): 3429-3452.  

[7] Cusworth, Daniel H., et al. "Using remote sensing to detect, validate, and quantify methane emissions from California solid waste operations." Environmental Research Letters 15.5 (2020): 054012.  

[8] Guha, Abhinav, et al. "Assessment of regional methane emission inventories through airborne quantification in the San Francisco Bay Area." Environmental Science & Technology 54.15 (2020): 9254-9264.  

[9] Krautwurst, Sven, et al. "Methane emissions from a Californian landfill, determined from airborne remote sensing and in situ measurements." Atmospheric Measurement Techniques 10.9 (2017): 3429-3452.  

[10] L. Guanter et al., "Mapping methane point emissions with the PRISMA spaceborne imaging spectrometer," Remote Sensing of Environment, vol. 265, p. 112671, 2021.  

[11] Nesme, Nicolas, et al. "Joint use of in-scene background radiance estimation and optimal estimation methods for quantifying methane emissions using PRISMA hyperspectral satellite data: Application to the Korpezhe industrial site." Remote Sensing 13.24 (2021): 4992.  

[12] C. C. Funk et al., "Clustering to improve matched filter detection of weak gas plumes in hyperspectral thermal imagery," IEEE Trans. Geosci. Remote Sens., vol. 39, no. 7, pp. 1410–1420, Jul. 2001, doi: 10.1109/36.934073.  

[13] Dowd E., Manning A. J., Orth-Lashley B. et al., First validation of high-resolution satellite-derived methane emissions from an active gas leak in the UK, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2023-2246, 2023.  

[14] Thorpe, A.K., Roberts, D.A., Bradley, E.S., Funk, C.C., Dennison, P.E., Leifer, I. (2013). "High resolution mapping of methane emissions from marine and terrestrial sources using a Cluster-Tuned Matched Filter technique and imaging spectrometry." Remote Sens. Environ., 134, 305–318.  

[15] Thompson, D.R., et al. "Real-time remote detection and measurement for airborne imaging spectroscopy: A case study with methane." Atmospheric Measurement Techniques 9.10 (2016): 5023-5033.  

[16] Niu, S., Golowich, S.E., Ingle, V.K., Manolakis, D.G. (2013). "New approach to remote gas-phase chemical quantification: Selected-band algorithm." Opt. Eng., 53, 021111.  

[17] Hulley, G.C., et al. "High spatial resolution imaging of methane and other trace gases with the Hyperspectral Thermal Emission Spectrometer (HyTES)." Atmospheric Measurement Techniques 9.5 (2016): 2393-2408.  

[18] Dennison, P.E., Thorpe, A.K., Roberts, D.A., Green, R.O. (2013). "Modeling sensitivity of imaging spectrometer data to carbon dioxide and methane plumes." In Proceedings of the Workshop on Hyperspectral Image and Signal Processing, Evolution in Remote Sensing, Gainesville, FL, USA, 25–28 June 2013; Volume 2013, pp. 46–49.  

[19] Theiler, J., Foy, B.R. (2006). Effect of Signal Contamination in Matched-Filter Detection of the Signal on a Cluttered Background. IEEE Geosci. Remote Sens. Lett., 3, 98–102.  
</small>
