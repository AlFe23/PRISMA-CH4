# PRISMA-CH4

CODE AVAILABLE FROM 7th JULY

This repository contains the code and resources for detecting and quantifying methane (CH4) emissions using PRISMA hyperspectral satellite imagery. 

The methodology, originally developed for CLEAR-UP project funded by ASI, leverages PRISMA Hyperspectral sensor to monitor methane emissions using an enhanced matched filter technique known as the Cluster Tuned Matched Filter (CTMF), applied to L1C radiance images. The unit absorption spectrum, obtained using multiple MODTRAN6® simulations, serves as the target in the matched filter, enabling methane detection and quantification. Key steps include utilizing PRISMA L1C radiance images in the 2300 nm short-wave infrared methane absorption window, adapting the matched filter to account for background clutter, and applying clustering for refined background spectrum calculations. Methane concentration is estimated by linearizing the Beer-Lambert absorption law and using the unit absorption spectrum, which represents the unit methane absorption as a function of wavelength. Automation through a lookup table (LUT) based on precomputed radiances for various conditions enables real-time applicability. Results computed over 200 PRISMA acquisitions revealed methane plumes at various sites of interest in countries including Algeria, Argentina, Brazil, China, India, Mexico, Turkmenistan, and the USA. For instance, the Natural Gas compression plant in Kamyshlydzha, Turkmenistan, consistently exhibited methane plumes across all 13 images analyzed. Several large landfills were examined using available PRISMA images or through specific acquisition requests. Evidence of methane concentration enhancement was observed in the Buenos Aires landfill and the Pirana landfill in Ahmedabad.

cite:
▪	Ferrari, A., Laneve, G., Pampanoni, V., Carvajal, A., Rossi, F. (2024, July). Monitoring Methane Emissions from Landfills Using PRISMA Imagery.  In IGARSS 2024-2024 IEEE International Geoscience and Remote Sensing Symposium (soon to be published). IEEE. (GitHub code Repository available here)



Links:

docker ctmf_v4_6:  https://drive.google.com/file/d/197rUulwsnqs67gWljOwPchdRw3-Yu87Y/view?usp=sharing


PRISMA sample imgs: https://drive.google.com/drive/folders/1IqJE_szLeWtHDRRdjURAl2JIBeQL_9zy?usp=sharing




METHODOLOGY
 
In the SWIR spectrum, the absorption windows of methane are centered about 1700nm and 2300nm. 
The radiance measured by a hyperspectral sensor is modeled as a superposition of signal b, scaled by its strength 𝛼, an average background radiance 𝑳_𝒎𝒆𝒂𝒏 (averaged over the entire image), and a zero-mean noise or clutter term ɛ [12][19].  
𝑟= 𝛼𝒃+𝑳_𝒎𝒆𝒂𝒏+ ɛ                (1)
 
ɛ is indicative of both sensor noise and scene clutter.( non-desirable components). For an uncorrelated background, the optimal filter is the basic matched filter (i.e. the target itself, scaled to make the variance of the filtered image equal to one).
Considering the reral case of a background clutter with correlation between spectral channels, the optimal filter is “matched” both to the signature of the target and the background: 
Clutter Matched Filter (CMF):     𝒒=(𝑪^(−1) 𝒃)/√(𝒃^𝑇 𝑪^(−1) 𝒃)                (2)
Here, 𝒒 is normalized to ensure that in the absence of signal, the variance of the matched filter image, 𝒒^𝑇 𝒓_𝒊 (for i-th pixel) is one. Values larger than one indicate strong evidence of the signature's presence, quantified as a "number of sigmas", Assuming that the clutter itself is Gaussian. When the clutter's covariance matrix 𝑪 is accurately known, equation (2) defines the optimal matched filter, maximizing the signal to clutter ratio. However, this covariance is seldom known beforehand and is typically deduced from the data, often by averaging the outer product of the mean-subtracted radiance across all pixels.
𝑪≈1/𝑁  ∑128_(𝑖=1)^𝑁▒〖(𝑳_𝒊−𝑳_𝒎𝒆𝒂𝒏 ) (𝑳_𝒊−𝑳_𝒎𝒆𝒂𝒏 )^𝑇          (3)〗

The CTMF [12], performs image clustering before applying the matched filter [14][15][16][17][18]. The average background spectrum and covariance matrix are calculated for each class. The CTMF score for the pixel (x,y) in the j-th class is given by:
𝛼_𝑗 (𝑥,𝑦)=𝒒_𝑗 (𝑳_((𝑥,𝑦)_𝑗 )−𝑳_(𝑚𝑒𝑎𝑛_𝑗 ) )   (4) 〖           𝒒〗_𝑗=(𝒃_𝑗^𝑇 𝑪_𝑗^(−1))/√(𝒃_𝑗^𝑇 𝑪_𝑗^(−1) 𝒃_𝑗 )            (5)
Where the optimal filter qj for the j-th class is expressed as [12],  𝑳_(𝑚𝑒𝑎𝑛,𝑗) is the average radiance of class j and 𝑪_𝒋 is the related covariance matrix [18] and 𝒃_𝒋 is the target signal.
 Where an important assumption is made approximating the off-plume radiance with the average radiance of class j, Lmean,j, to determine the expression of b.  

𝑳_(𝑛𝑜_𝑝𝑙)=𝑳_(𝑚𝑒𝑎𝑛,𝑗)     (8)
![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/0bac7d17-1504-4945-8854-5a8a6ebce566)



METHANE CONCENTRATION ESTIMATION	
 
The impact of methane enhancement is described by the Beer-Lambert absorption law; linearization of Beer-Lambert law using a first-order Taylor series approximation allows for the formulation of a linearion inversion problem to estimate methane concentration:
 
𝜌_𝐿𝑀=((𝑳_((𝑥,𝑦)_𝑗 )−𝑳_(𝑚𝑒𝑎𝑛_𝑗 ) ) ⅆ_𝑗^𝑇 𝑪_𝑗^(−1))/(𝒅_𝑗^𝑇 𝑪_𝑗^(−1) 𝒅_𝑗 )          (8) 〖         ⅆ〗_𝑗=−𝐴(1+1/cos⁡𝜗 ) 𝑳_(𝑚𝑒𝑎𝑛,𝑗)	
(9) 
Where θ is the solar zenith angle, A is the absorption coefficient at the instrument's wavelengths, representing the specific absorption of methane.
The linear estimation method allows to calculate the concentration of the gas based on the difference between the observed radiance 𝑳_((𝑥_𝑖,𝑦_𝑖 ) ) and the mean radiance 𝑳_(𝑚𝑒𝑎𝑛,𝑗). This approach is consistent with the methodology outlined in previous studies [11] and allows to calculate the concentration of the gas as the CTMF score 𝛼_𝑗 divided by the optimal filter applied to the target signal (𝒒_𝑗∙𝒃_𝑗).
![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/cc88832b-740c-463a-b6f9-3c4b71f7340c)



Scene Specific Unit Target Spectrum
 
The unit absorption spectrum is derived from multiple simulations with varying methane gas concentrations using radiative transfer models such as MODTRAN6®.  This spectrum represents the methane absorption as a function of wavelength, normalized per unit path length and concentration. The unit for this measurement is typically (ppm·m)-1, indicating the inverse of the product of methane concentration (in parts per million) and the path length (in meters) over which the absorption occurs.
Subsequently the simulated radiance spectra are convolved with PRISMA spectral response function, and a regression analysis is carried out to obtain the unit absorption spectrum for each PRISMA spectral channel. The slope of the regression line, which best represents the relationship between the concentration-path product and the natural logarithm of the radiance, provides the unit absorption value for each specific wavelength band.  This unit absorption spectrum is scaled by element-wise multiplication with the average radiance at each wavelength to generate the target spectrum, which is then used in the matched filter. To automate the generation of the unit target spectrum and make it independent of real-time MODTRAN6® runs, a lookup table (LUT) has been developed based on the methodology described by Foote et al. (2020) [4]. This LUT contains precomputed at-sensor radiances for a wide range of atmospheric and geometric parameters, including variations in sensor altitude, water vapor content, ground elevation, sun zenith angle, methane concentration. For each PRISMA acquisition, the Level 2C PRISMA image is also used to read the water vapor content, and a Digital Elevation Model (DEM) is used to read the ground elevation. By interpolating within the LUT, it is possible to quickly generate a unit target spectrum that accurately matches the conditions of a given satellite pass.
![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/85df8a51-21d2-403e-9e61-2eeb3e1846d3)




![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/631d12fe-f2eb-424c-9e27-da3ace35469a)


![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/160778eb-f03e-477b-83be-0ce0200c0bbb)


![image](https://github.com/AlFe23/PRISMA-CH4/assets/105355911/d471e0b0-3343-4c26-a5ef-7cb7f16bfa4d)


