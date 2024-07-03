# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 11:00:29 2024

@author: Alvise Ferrari

@description: 
    The methodology, originally developed for CLEAR-UP project funded by ASI, leverages PRISMA Hyperspectral sensor to monitor methane 
    emissions using an enhanced matched filter technique known as the Cluster Tuned Matched Filter (CTMF), applied to L1C radiance images. 
    The unit absorption spectrum, obtained using multiple MODTRAN6® simulations, serves as the target in the matched filter, enabling methane 
    detection and quantification. Key steps include utilizing PRISMA L1C radiance images in the 2300 nm short-wave infrared methane absorption 
    window, adapting the matched filter to account for background clutter, and applying clustering for refined background spectrum calculations. 
    Methane concentration is estimated by linearizing the Beer-Lambert absorption law and using the unit absorption spectrum, which represents 
    the unit methane absorption as a function of wavelength. Automation through a lookup table (LUT) based on precomputed radiances for various 
    conditions enables real-time applicability.
"""





import numpy as np
import os
import matplotlib.pyplot as plt
import h5py
from sklearn.cluster import KMeans
from numpy.linalg import inv
from osgeo import gdal, osr
import sys
from os.path import exists
import scipy.ndimage
import argparse
import spectral
import time
import xarray as xr
import pandas as pd
import argparse


import numpy as np
from osgeo import gdal
import subprocess
import os


def prismaL2C_WV_read(filename):
    
    # Open the HDF5 file
    with h5py.File(filename, 'r') as f:
        # Read the Water vapor Mask from L2C PRISMA data
        PRS_L2C_WVM = f['HDFEOS/SWATHS/PRS_L2C_WVM/Data Fields/WVM_Map'][:]
        latitude_WVM = f['HDFEOS/SWATHS/PRS_L2C_WVM/Geolocation Fields/Latitude']
        longitude_WVM = f['HDFEOS/SWATHS/PRS_L2C_WVM/Geolocation Fields/Longitude']
        
        # Read scale factors to transform DN (uint16) into WV physical values (g/cm2)
        L2ScaleWVMMin = f.attrs['L2ScaleWVMMin']
        L2ScaleWVMMax = f.attrs['L2ScaleWVMMax']
       
    #Convert DN to WV in g/cm^2    (pg.213 product spec doc)
    WVM_unit = L2ScaleWVMMin + PRS_L2C_WVM*(L2ScaleWVMMax-L2ScaleWVMMin)/65535
    meanWV = np.mean(WVM_unit)
    # Print the mean Water Vapor value
    print("Mean Water Vapor (g/cm^2):", meanWV)
    
    return  meanWV, PRS_L2C_WVM, latitude_WVM, longitude_WVM

def prismaL1_SZA_read(filename):
    
    # Open the HDF5 file
    with h5py.File(filename, 'r') as f:

        # Read SZA in degrees
        SZA = f.attrs['Sun_zenith_angle']

    # Print the Sun Zenith Angle
    print("Sun Zenith Angle (degrees):", SZA)
    
    return  SZA

def prismaL2C_bbox_read(filename):
    # Open the HDF5 file
    with h5py.File(filename, 'r') as f:
        # Read the geolocation fields
        latitude_WVM = f['HDFEOS/SWATHS/PRS_L2C_WVM/Geolocation Fields/Latitude'][:]
        longitude_WVM = f['HDFEOS/SWATHS/PRS_L2C_WVM/Geolocation Fields/Longitude'][:]

    # Calculate bounding box
    min_lat = np.min(latitude_WVM)
    max_lat = np.max(latitude_WVM)
    min_lon = np.min(longitude_WVM)
    max_lon = np.max(longitude_WVM)
    
    return (min_lon, max_lon, min_lat, max_lat)

def mean_elev_fromDEM(dem_file, bbox):
    # Load the NetCDF file
    ds = xr.open_dataset(dem_file)
    
    # Extract bounding box coordinates
    min_lon, max_lon, min_lat, max_lat = bbox
    
    # Slice the dataset to only include the area within the bounding box
    elevation_subset = ds.sel(lon=slice(min_lon, max_lon), lat=slice(min_lat, max_lat))
    
    # Calculate and print mean elevation within the bounding box expressed in Km
    mean_elevation_subset = elevation_subset['elev'].mean() / 1000
    print("Mean Elevation within Bounding Box in Km:", mean_elevation_subset.values)

    # # Convert elevation data to a DataFrame (optional)
    # elevation_df = elevation_subset['elev'].to_dataframe()
    # print(elevation_df.head())

    # # Plot elevation subset
    # elevation_subset['elev'].plot()
    # plt.show()

    # Close the dataset
    ds.close()
    return mean_elevation_subset.values

def prisma_read(filename):
    
    # Open the HDF5 file
    with h5py.File(filename, 'r') as f:
        # Read the VNIR and SWIR spectral data
        vnir_cube_DN = f['HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/VNIR_Cube'][:]
        swir_cube_DN = f['HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SWIR_Cube'][:]
    
        # Read central wavelengths and FWHM for VNIR and SWIR bands
        cw_vnir = f['KDP_AUX/Cw_Vnir_Matrix'][:]
        fwhm_vnir = f['KDP_AUX/Fwhm_Vnir_Matrix'][:]
        cw_swir = f['KDP_AUX/Cw_Swir_Matrix'][:]
        fwhm_swir = f['KDP_AUX/Fwhm_Swir_Matrix'][:]
    
        # Read scale factor and offset to transform DN to radiance physical value
        offset_swir = f.attrs['Offset_Swir']
        scaleFactor_swir = f.attrs['ScaleFactor_Swir']
        offset_vnir = f.attrs['Offset_Vnir']
        scaleFactor_vnir = f.attrs['ScaleFactor_Vnir']
        
        # Read geolocation fields
        latitude_vnir = f['HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Latitude_VNIR'][:]
        longitude_vnir = f['HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Longitude_VNIR'][:]
        latitude_swir = f['HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Latitude_SWIR'][:]
        longitude_swir = f['HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Longitude_SWIR'][:]
    
    # Convert DN to radiance for VNIR and SWIR data #[W/(str*um*m^2)]
    swir_cube_rads = (swir_cube_DN / scaleFactor_swir) - offset_swir
    vnir_cube_rads = (vnir_cube_DN / scaleFactor_vnir) - offset_vnir
    
    #Convert PRISMA radiance from [W/(str*um*m^2)] to [μW*cm-2*nm-1*sr-1] in order to meet original AVIRIS radiance unit.
    swir_cube_rads = swir_cube_rads * 0.1
    vnir_cube_rads = vnir_cube_rads * 0.1
    
    # The variables swir_cube and vnir_cube now contain the physical radiance values
    
    
    #############################################################################
    # Extract the specific bands for RGB
    red_rad_channel = vnir_cube_rads[:, 29, :]   # Red channel
    green_rad_channel = vnir_cube_rads[:, 19, :] # Green channel
    blue_rad_channel = vnir_cube_rads[:, 7, :]   # Blue channel
    
    # Normalize each channel to the range [0, 1]
    red_rad_normalized = (red_rad_channel - red_rad_channel.min()) / (red_rad_channel.max() - red_rad_channel.min())
    green_rad_normalized = (green_rad_channel - green_rad_channel.min()) / (green_rad_channel.max() - green_rad_channel.min())
    blue_rad_normalized = (blue_rad_channel - blue_rad_channel.min()) / (blue_rad_channel.max() - blue_rad_channel.min())
    
    # Combine the channels into an RGB image
    rgb_image = np.stack([red_rad_normalized, green_rad_normalized, blue_rad_normalized], axis=-1)
    
    # # Create a larger figure
    # plt.figure(figsize=(10, 10))  # You can adjust the size as needed
    
    # # Plotting
    # plt.imshow(rgb_image)
    # plt.title('RGB Image from Radiance Data')
    # plt.axis('off')  # Turn off axis numbers and labels
    # plt.show(block=False)
    # plt.show()
    #############################################################################
    
    # Convert from BIL to BIP format ( BIL = M-3-N ; BIP = 3-M-N ; BSQ = M-N-3 )
    vnir_cube_bip = np.transpose(vnir_cube_rads, (0, 2, 1))
    swir_cube_bip = np.transpose(swir_cube_rads, (0, 2, 1))
    
    # Rotate 270 degrees counterclockwise (equivalent to 90 degrees clockwise) and flip horizontally
    vnir_cube_bip = np.rot90(vnir_cube_bip, k=-1, axes=(0, 1))  # Rotate 270 degrees counterclockwise
    #vnir_cube_bip = np.flip(vnir_cube_bip, axis=1)  # Flip horizontally along the columns axis
    swir_cube_bip = np.rot90(swir_cube_bip, k=-1, axes=(0, 1))  # Rotate 270 degrees counterclockwise
    #swir_cube_bip = np.flip(swir_cube_bip, axis=1)  # Flip horizontally along the columns axis
    
    # The variables swir_cube_bip and vnir_cube_bip now contain the physical radiance values in BIP format
    
    # PROBLEMA:
    # Si è trovato che le bande 1,2,3 del VNIR cube e le bande 172,173 dello SWIR cube sono azzerate:
    # Q&A PRISMA_ATBD.pdf pg.118	
    # VNIR: Remove bands 1, 2, 3 (0-indexed: 0, 1, 2), for these bands CWs and FWHMs are already set to 0.
    VNIR_cube_clean = np.delete(vnir_cube_bip, [0, 1, 2], axis=2)
    # SWIR: Remove bands 172, 173 (0-indexed: 171, 172), for these bands CWs and FWHMs are already set to 0.
    SWIR_cube_clean = np.delete(swir_cube_bip, [171, 172], axis=2)
    # SWIR: Remove bands 1, 2, 3, 4, since they are redundant with the last 4 bands of VNIR cube
    SWIR_cube_clean = np.delete(SWIR_cube_clean, [0, 1, 2, 3], axis=2)
    
    #Extract actual values of CWs and FWHMs. They are stored in standars (1000,256) arrays and have to be extracted
    cw_vnir = cw_vnir[:, 99:162]
    fwhm_vnir = fwhm_vnir[:, 99:162]
    cw_swir = cw_swir[:, 81:252]
    fwhm_swir = fwhm_swir[:, 81:252]
    
    #Reverse CWs and FWHMs vectors as they are in decrasing frequency order, but we want them opposite
    cw_vnir = cw_vnir[:, ::-1]
    fwhm_vnir = fwhm_vnir[:, ::-1]
    cw_swir = cw_swir[:, ::-1]
    fwhm_swir = fwhm_swir[:, ::-1]
    
    # SWIR: Remove bands 1, 2, 3, 4, since they are redundant with the last 4 bands of VNIR cube
    cw_swir_clean = np.delete(cw_swir, [0, 1, 2, 3], axis=1)
    fwhm_swir_clean = np.delete(fwhm_swir, [0, 1, 2, 3], axis=1)
    
    
    #Let's now concatenate the arrays for radiance cube, central wavelengths and FWHMs
    # Concatenate VNIR and SWIR cubes along the band direction
    # Make sure VNIR_cube_filtered and SWIR_cube_filtered are your actual data cubes
    concatenated_cube = np.concatenate((SWIR_cube_clean, VNIR_cube_clean), axis=2)
    # Reverse frequnecies representation according to CWs and FWHMs vectors
    concatenated_cube = concatenated_cube[:, :, ::-1]
    #####################################concatenated_cube =  np.rot90(concatenated_cube,k=3) # si è spostata questa operazione ai due cubi VNIR e SWIR direttamente
    # Concatenate CW and FWHM arrays
    # These should be the actual SRF variables
    concatenated_cw = np.concatenate((cw_vnir, cw_swir_clean), axis=1)
    concatenated_fwhm = np.concatenate((fwhm_vnir, fwhm_swir_clean), axis=1)
    
    #############################################################################
    #Plotting radiance at each band, averaged over all pixels of the image
    #let's plot the average radiance
    mean_rads_concatenated_cube = np.mean(concatenated_cube, axis=(0,1))
    
    # # Coordinates of the pixel (replace with your desired coordinates)
    # x, y = 700, 300
    # # Extracting the radiance values for the selected pixel
    # radiance_values = concatenated_cube[x, y, :]
    
    # # Plotting 
    # plt.figure(figsize=(10, 6))
    # plt.plot(np.mean(concatenated_cw, axis=0), mean_rads_concatenated_cube, label='Radiance vs Wavelength')
    # plt.xlabel('Central Wavelength (nm)')
    # plt.ylabel('Radiance')
    # plt.title('average Radiance Spectrum (vs. mean CWs)')
    # plt.legend()
    # plt.show(block=False)
    # plt.show()
    #############################################################################
    
    #Print the RGB image
    # Assuming rads_array is your (1000, 1000, n_bands) array with radiance data

    # Extract the specific bands for RGB
    red_channel = concatenated_cube[:, :, 29]   # Red channel
    green_channel = concatenated_cube[:, :, 19] # Green channel
    blue_channel = concatenated_cube[:, :, 7]   # Blue channel
    
    # Normalize each channel to the range [0, 1]
    red_normalized = (red_channel - red_channel.min()) / (red_channel.max() - red_channel.min())
    green_normalized = (green_channel - green_channel.min()) / (green_channel.max() - green_channel.min())
    blue_normalized = (blue_channel - blue_channel.min()) / (blue_channel.max() - blue_channel.min())
    
    # Combine the channels into an RGB image
    rgb_image = np.stack([red_normalized, green_normalized, blue_normalized], axis=-1)
    
    # # Create a larger figure
    # plt.figure(figsize=(10, 10))  # You can adjust the size as needed
    
    # # Plotting
    # plt.imshow(rgb_image)
    # plt.title('RGB Image from Radiance Data')
    # plt.axis('off')  # Turn off axis numbers and labels
    # plt.show(block=False)
    # plt.show()
    #############################################################################
    
    return concatenated_cube, concatenated_cw, concatenated_fwhm, rgb_image, vnir_cube_bip, swir_cube_bip, latitude_vnir, longitude_vnir, latitude_swir, longitude_swir

def k_means_hyperspectral(image, k):
    """
    Parameters:
    image (numpy array): L'immagine iperspettrale con dimensioni (n_righe, n_colonne, n_bande).
    k (int): Il numero di cluster da utilizzare nel K-means.

    Returns:
    numpy array: Immagine classificata con dimensioni (n_righe, n_colonne).
    """

    # Ottieni le dimensioni dell'immagine
    n_righe, n_colonne, n_bande = image.shape

    # Rimodella l'immagine per il clustering
    reshaped_image = image.reshape((n_righe * n_colonne, n_bande))

    # Applica il K-means clustering
    kmeans = KMeans(n_clusters=k, n_init=10)
    kmeans.fit(reshaped_image)

    # Ottieni le etichette dei cluster e rimodella per l'immagine classificata
    clustered_image = kmeans.labels_.reshape((n_righe, n_colonne))

    return clustered_image

def plot_classified_image(classified_image, title='Classified Image'):
    """
    Visualizza l'immagine classificata.

    Parametri:
    classified_image (numpy array): Immagine classificata con dimensioni (n_righe, n_colonne).
    title (str): Titolo del grafico.
    """

    plt.figure(figsize=(8, 6))
    plt.imshow(classified_image, cmap='jet')
    plt.colorbar()
    plt.title(title)
    plt.show(block=False)
    plt.show()

def calculate_statistics(image, classified_image, k):
    """
    Calcola la radianza media e la matrice di covarianza per ogni cluster.

    Parametri:
    image (numpy array): Immagine iperspettrale originale con dimensioni (n_righe, n_colonne, n_bande).
    classified_image (numpy array): Immagine classificata con dimensioni (n_righe, n_colonne).
    k (int): Numero di cluster.

    Ritorna:
    mean_radiance (list of numpy arrays): Radianza media per ogni cluster.
    covariance_matrices (list of numpy arrays): Matrici di covarianza per ogni cluster.
    """

    n_righe, n_colonne, n_bande = image.shape
    mean_radiance = []
    covariance_matrices = []

    # Rimodella l'immagine per un facile accesso ai pixel
    reshaped_image = image.reshape((n_righe * n_colonne, n_bande))

    for cluster in range(k):
        # Estrai i pixel di questo cluster
        pixels = reshaped_image[classified_image.reshape(n_righe * n_colonne) == cluster]

        # Calcola la radianza media
        mean_radiance.append(np.mean(pixels, axis=0))

        # Calcola la matrice di covarianza
        covariance_matrices.append(np.cov(pixels, rowvar=False))

    return mean_radiance, covariance_matrices

def extract_number(file_path):
    file_name = os.path.basename(file_path)  # Extract the file name from the full path
    return int(file_name.split('.')[0])

def generate_template_from_bands(centers, fwhm, simRads, simWave, concentrations):

    """Calculate a unit absorption spectrum for methane by convolving with given band information.

    :param centers: wavelength values for the band centers, provided in nanometers.
    :param fwhm: full width half maximum for the gaussian kernel of each band.
    :return template: the unit absorption spectum
    """
    # import scipy.stats
    #SCALING = 1e5
    SCALING = 1
    # centers = np.asarray(centers)
    # fwhm = np.asarray(fwhm)
    if np.any(~np.isfinite(centers)) or np.any(~np.isfinite(fwhm)):
        raise RuntimeError('Band Wavelengths Centers/FWHM data contains non-finite data (NaN or Inf).')
    if centers.shape[0] != fwhm.shape[0]:
        raise RuntimeError('Length of band center wavelengths and band fwhm arrays must be equal.')
    # lib = spectral.io.envi.open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ch4.hdr'),
    #                             os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ch4.lut'))
    # ######## Run these 3 line instead of the command above if you want to run section of the code with 'F9'
    # current_dir = os.getcwd()
    # lib = spectral.io.envi.open(os.path.join(current_dir, 'ch4.hdr'),
    #                             os.path.join(current_dir, 'ch4.lut'))
    # ####################################################################
    # rads = np.asarray(lib.asarray()).squeeze()
    # wave = np.asarray(lib.bands.centers)
    #concentrations = np.asarray([0, 500, 1000, 2000, 4000, 8000, 16000])
    # sigma = fwhm / ( 2 * sqrt( 2 * ln(2) ) )  ~=  fwhm / 2.355
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    # response = scipy.stats.norm.pdf(wave[:, None], loc=centers[None, :], scale=sigma[None, :])
    # Evaluate normal distribution explicitly
    var = sigma ** 2
    denom = (2 * np.pi * var) ** 0.5
    numer = np.exp(-(simWave[:, None] - centers[None, :])**2 / (2*var))
    response = numer / denom
    # Normalize each gaussian response to sum to 1.
    response = np.divide(response, response.sum(axis=0), where=response.sum(axis=0) > 0, out=response)
    # implement resampling as matrix multiply
    resampled = simRads.dot(response)
    lograd = np.log(resampled, out=np.zeros_like(resampled), where=resampled > 0)
    slope, _, _, _ = np.linalg.lstsq(np.stack((np.ones_like(concentrations), concentrations)).T, lograd, rcond=None)
    spectrum = slope[1, :] * SCALING
    target = np.stack((centers, spectrum)).T  # np.stack((np.arange(spectrum.shape[0]), centers, spectrum)).T
    return target

def calculate_matched_filter(rads_array, classified_image, mean_radiance, covariance_matrices, target_spectra, k):
   
    n_rows, n_columns, n_bands = rads_array.shape
    matched_filter_scores = np.zeros((n_rows, n_columns))
    concentration_map = np.zeros((n_rows, n_columns))
    
    for w in range(n_columns):
        # Estrarre lo spettro target per la w-esima colonna
        target_spectrum = target_spectra[:, w]

        # Calcolare lo spettro target specifico per ogni classe alla w-esima colonna
        class_target_spectra = [target_spectrum * mean_radiance[cls] for cls in range(k)]

        # Calcolare il filtro adattato per ogni classe alla w-esima colonna
        adapted_filters = []
        for cls in range(k):
            target_spc_cls = class_target_spectra[cls]
            inv_cov_matrix_cls = inv(covariance_matrices[cls])
            q_cls = np.dot(target_spc_cls.T, inv_cov_matrix_cls)
            q_cls /= np.sqrt(np.dot(np.dot(target_spc_cls.T, inv_cov_matrix_cls), target_spc_cls)) # il termine al denominatore dä un numero, quindi si puö fare divisione semplice
            adapted_filters.append(q_cls)

        for i in range(n_rows):
            # Determinare la classe del pixel
            pixel_class = classified_image[i, w]

            # Calcolare lo score del matched filter per il pixel
            pixel_spectrum = rads_array[i, w]
            matched_filter_score = np.dot(adapted_filters[pixel_class], pixel_spectrum - mean_radiance[pixel_class])
            matched_filter_scores[i, w] = matched_filter_score
            concentration_map[i,w] = matched_filter_score/np.dot(q_cls,target_spc_cls)

    return matched_filter_scores, concentration_map



def plot_matched_filter_scores(scores, title='Matched Filter Scores'):
    """
    Visualizza i punteggi del filtro adattato utilizzando la palette turbo.

    Parametri:
    scores (numpy array): Array dei punteggi del filtro adattato.
    title (str): Titolo del grafico.
    """

    plt.figure(figsize=(10, 10))
    plt.imshow(scores, cmap='turbo')
    plt.colorbar()
    plt.title(title)
    plt.show(block=False)
    plt.show()

def save_as_geotiff_rgb(rgb_data, output_file, latitude_vnir, longitude_vnir):
    """
    Save a numpy RGB array as a GeoTIFF file with proper georeferencing using gdalwarp for accurate projection.
    
    Parameters:
    rgb_data (numpy array): RGB data to save.
    output_file (str): Path to the output file.
    latitude_vnir (numpy array): Array of latitude values.
    longitude_vnir (numpy array): Array of longitude values.
    """
    # Ensure the data is in float32 format
    rgb_data_float32 = rgb_data.astype(np.float32)

    # Revert the rotation applied in prisma_read
    rgb_data_float32 = np.rot90(rgb_data_float32, k=1, axes=(0, 1))

    # Temporary files to hold the intermediate GeoTIFF and VRT
    temp_file = 'temp_output.tif'
    vrt_file = 'temp_output.vrt'
    lat_file = 'latitude.tif'
    lon_file = 'longitude.tif'

    # Create temporary files for latitude and longitude
    driver = gdal.GetDriverByName('GTiff')
    lat_ds = driver.Create(lat_file, latitude_vnir.shape[1], latitude_vnir.shape[0], 1, gdal.GDT_Float32)
    lon_ds = driver.Create(lon_file, longitude_vnir.shape[1], longitude_vnir.shape[0], 1, gdal.GDT_Float32)

    lat_ds.GetRasterBand(1).WriteArray(latitude_vnir.astype(np.float32))
    lon_ds.GetRasterBand(1).WriteArray(longitude_vnir.astype(np.float32))

    lat_ds = None
    lon_ds = None

    # Create a temporary GeoTIFF file without detailed geotransform and projection
    dataset = driver.Create(temp_file, rgb_data.shape[1], rgb_data.shape[0], 3, gdal.GDT_Float32)
    
    # Write the RGB data for each band
    for i in range(3):
        dataset.GetRasterBand(i + 1).WriteArray(rgb_data_float32[:, :, i])

    # Save and close the dataset
    dataset.FlushCache()
    dataset = None

    # Create the VRT file that uses the latitude and longitude arrays for georeferencing
    vrt_options = gdal.TranslateOptions(format='VRT')
    gdal.Translate(vrt_file, temp_file, options=vrt_options)

    # Open the VRT file and set geolocation metadata
    vrt_ds = gdal.Open(vrt_file, gdal.GA_Update)
    vrt_ds.SetMetadata({
        'X_DATASET': lon_file,
        'X_BAND': '1',
        'Y_DATASET': lat_file,
        'Y_BAND': '1',
        'PIXEL_OFFSET': '0',
        'LINE_OFFSET': '0',
        'PIXEL_STEP': '1',
        'LINE_STEP': '1'
    }, 'GEOLOCATION')

    # Metadata statement
    description = ("This product has been generated by Alvise Ferrari for School of Aerospace Engineering, "
                          "La Sapienza, under terms of license of CLEAR-UP, a project funded by the Italian Space Agency. "
                          "the dissemination of this product is closely linked to the agreements established under the CLEAR-UP project. "
                          "The authors of the code by which the product was generated cannot be held responsible for any improper use or dissemination of this product")
    
    # Set specific metadata
    vrt_ds.SetMetadataItem('DESCRIPTION', description)
    vrt_ds = None

    # Use gdalwarp to apply accurate projection and georeferencing
    subprocess.run([
        'gdalwarp',
        '-geoloc',
        vrt_file,
        output_file
    ], check=True)

    # Remove the temporary files
    os.remove(temp_file)
    os.remove(vrt_file)
    os.remove(lat_file)
    os.remove(lon_file)

def save_as_geotiff_single_band(data, output_file, latitude_vnir, longitude_vnir):
    """
    Save a numpy single-band array as a GeoTIFF file with proper georeferencing using gdalwarp for accurate projection.

    Parameters:
    data (numpy array): Data to save.
    output_file (str): Path to the output file.
    latitude_vnir (numpy array): Array of latitude values.
    longitude_vnir (numpy array): Array of longitude values.
    """
    # Ensure the data is in float32 format
    data_float32 = data.astype(np.float32)

    # Revert the rotation applied in prisma_read
    data_float32 = np.rot90(data_float32, k=1, axes=(0, 1))

    # Temporary files to hold the intermediate GeoTIFF and VRT
    temp_file = 'temp_output.tif'
    vrt_file = 'temp_output.vrt'
    lat_file = 'latitude.tif'
    lon_file = 'longitude.tif'

    # Create temporary files for latitude and longitude
    driver = gdal.GetDriverByName('GTiff')
    lat_ds = driver.Create(lat_file, latitude_vnir.shape[1], latitude_vnir.shape[0], 1, gdal.GDT_Float32)
    lon_ds = driver.Create(lon_file, longitude_vnir.shape[1], longitude_vnir.shape[0], 1, gdal.GDT_Float32)

    lat_ds.GetRasterBand(1).WriteArray(latitude_vnir.astype(np.float32))
    lon_ds.GetRasterBand(1).WriteArray(longitude_vnir.astype(np.float32))

    lat_ds = None
    lon_ds = None

    # Create a temporary GeoTIFF file without detailed geotransform and projection
    dataset = driver.Create(temp_file, data.shape[1], data.shape[0], 1, gdal.GDT_Float32)
    
    # Write the data
    dataset.GetRasterBand(1).WriteArray(data_float32)

    # Save and close the dataset
    dataset.FlushCache()
    dataset = None

    # Create the VRT file that uses the latitude and longitude arrays for georeferencing
    vrt_options = gdal.TranslateOptions(format='VRT')
    gdal.Translate(vrt_file, temp_file, options=vrt_options)

    # Open the VRT file and set geolocation metadata
    vrt_ds = gdal.Open(vrt_file, gdal.GA_Update)
    vrt_ds.SetMetadata({
        'X_DATASET': lon_file,
        'X_BAND': '1',
        'Y_DATASET': lat_file,
        'Y_BAND': '1',
        'PIXEL_OFFSET': '0',
        'LINE_OFFSET': '0',
        'PIXEL_STEP': '1',
        'LINE_STEP': '1'
    }, 'GEOLOCATION')

    # Metadata statement
    description = ("This product has been generated by Alvise Ferrari for School of Aerospace Engineering, "
                          "La Sapienza, under terms of license of CLEAR-UP, a project funded by the Italian Space Agency. "
                          "the dissemination of this product is closely linked to the agreements established under the CLEAR-UP project. "
                          "The authors of the code by which the product was generated cannot be held responsible for any improper use or dissemination of this product")
    
    # Set specific metadata
    vrt_ds.SetMetadataItem('DESCRIPTION', description)
    vrt_ds = None

    # Use gdalwarp to apply accurate projection and georeferencing
    subprocess.run([
        'gdalwarp',
        '-geoloc',
        vrt_file,
        output_file
    ], check=True)

    # Remove the temporary files
    os.remove(temp_file)
    os.remove(vrt_file)
    os.remove(lat_file)
    os.remove(lon_file)
    

###################################################################################
#FUNZIONI PER ESTRARRE E DEFINIRE LE VARIABILI NECESSARIE ALL'ESTRAZIONE DELLE FIRME SPETTRALI DALLA LUT



def check_param(value, min, max, name):
    if value < min or value > max:
        raise ValueError('The value for {name} exceeds the sampled parameter space.'
                         'The limits are[{min}, {max}], requested {value}.')

@np.vectorize
# [0.,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80]
def get_5deg_zenith_angle_index(zenith_value):
    check_param(zenith_value, 0, 80, 'Zenith Angle')
    return zenith_value / 5

@np.vectorize
def get_5deg_sensor_height_index(sensor_value):  # [1, 2, 4, 10, 20, 120]
    # Only check lower bound here, atmosphere ends at 120 km so clamping there is okay.
    check_param(sensor_value, 1, np.inf, 'Sensor Height')
    # There's not really a pattern here, so just linearly interpolate between values -- piecewise linear
    if sensor_value < 1.0:
        return np.float64(0.0)
    elif sensor_value < 2.0:
        idx = sensor_value - 1.0
        return idx
    elif sensor_value < 4:
        return sensor_value / 2
    elif sensor_value < 10:
        return (sensor_value / 6) + (4.0 / 3.0)
    elif sensor_value < 20:
        return (sensor_value / 10) + 2
    elif sensor_value < 120:
        return (sensor_value / 100) + 3.8
    else:
        return 5

@np.vectorize
def get_5deg_ground_altitude_index(ground_value):  # [0, 0.5, 1.0, 2.0, 3.0]
    check_param(ground_value, 0, 3, 'Ground Altitude')
    if ground_value < 1:
        return 2 * ground_value
    else:
        return 1 + ground_value

@np.vectorize
def get_5deg_water_vapor_index(water_value):  # [0,1,2,3,4,5,6]
    check_param(water_value, 0, 6, 'Water Vapor')
    return water_value

@np.vectorize
# [0.0,1000,2000,4000,8000,16000,32000,64000]
def get_5deg_methane_index(methane_value):
    # the parameter clamps should rarely be calle because there are default concentrations, but the --concentraitons parameter exposes these
    check_param(methane_value, 0, 64000, 'Methane Concentration')
    if methane_value <= 0:
        return 0
    elif methane_value < 1000:
        return methane_value / 1000
    return np.log2(methane_value / 500)

@np.vectorize
def get_carbon_dioxide_index(coo_value):
    check_param(coo_value, 0, 1280000, 'Carbon Dioxode Concentration')
    if coo_value <= 0:
        return 0
    elif coo_value < 20000:
        return coo_value / 20000
    return np.log2(coo_value / 10000)

def get_5deg_lookup_index(zenith=0, sensor=120, ground=0, water=0, conc=0, gas='ch4'):
    if 'ch4' in gas:
        idx = np.asarray([[get_5deg_zenith_angle_index(zenith)],
                          [get_5deg_sensor_height_index(sensor)],
                          [get_5deg_ground_altitude_index(ground)],
                          [get_5deg_water_vapor_index(water)],
                          [get_5deg_methane_index(conc)]])
    elif 'co2' in gas:
        idx = np.asarray([[get_5deg_zenith_angle_index(zenith)],
                          [get_5deg_sensor_height_index(sensor)],
                          [get_5deg_ground_altitude_index(ground)],
                          [get_5deg_water_vapor_index(water)],
                          [get_carbon_dioxide_index(conc)]])
    else:
        raise ValueError('Unknown gas provided.')
    return idx

def spline_5deg_lookup(grid_data, zenith=0, sensor=120, ground=0, water=0, conc=0, gas='ch4', order=1):
    coords = get_5deg_lookup_index(
        zenith=zenith, sensor=sensor, ground=ground, water=water, conc=conc, gas=gas)
    # correct_lookup = np.asarray([scipy.ndimage.map_coordinates(
    #     im, coordinates=coords, order=order, mode='nearest') for im in np.moveaxis(grid_data, 5, 0)])
    if order == 1:
        coords_fractional_part, coords_whole_part = np.modf(coords)
        #coords_near_slice = tuple((slice(int(c), int(c+2)) for c in coords_whole_part))  #This line gives a warning, so it as been modified as below
        coords_near_slice = tuple((slice(int(c[0]), int(c[0]+2)) if isinstance(c, np.ndarray) else slice(int(c), int(c+2)) for c in coords_whole_part))
        near_grid_data = grid_data[coords_near_slice]
        new_coord = np.concatenate((coords_fractional_part * np.ones((1, near_grid_data.shape[-1])),
                                    np.arange(near_grid_data.shape[-1])[None, :]), axis=0)
        lookup = scipy.ndimage.map_coordinates(near_grid_data, coordinates=new_coord, order=1, mode='nearest')
    elif order == 3:
        lookup = np.asarray([scipy.ndimage.map_coordinates(
            im, coordinates=coords_fractional_part, order=order, mode='nearest') for im in np.moveaxis(near_grid_data, 5, 0)])
    return lookup.squeeze()


###################################################################################
#FUNZIONI PER LA LETTURA DELLA LUT

def load_ch4_dataset(lut_file_path):
    # Ensure the function uses the passed file path instead of a hardcoded one
    datafile = h5py.File(lut_file_path, 'r', rdcc_nbytes=4194304)
    return datafile['modtran_data'], datafile['modtran_param'], datafile['wave'], 'ch4'

def generate_library(gas_concentration_vals, lut_file, zenith=0, sensor=120, ground=0, water=0, order=1, dataset_fcn=load_ch4_dataset):
    # Use the passed `dataset_fcn` function, allowing for flexibility in data loading.
    grid, params, wave, gas = dataset_fcn(lut_file)
    rads = np.empty((len(gas_concentration_vals), grid.shape[-1]))
    for i, ppmm in enumerate(gas_concentration_vals):
        rads[i, :] = spline_5deg_lookup(
            grid, zenith=zenith, sensor=sensor, ground=ground, water=water, conc=ppmm, gas=gas, order=order)
    return rads, np.array(wave)


###################################################################################



# Define the main processing function
def ch4_detection(L1_file, L2C_file, dem_file, lut_file, output_dir, k=30):
   
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
   # Read SZA from L1 image attributes
    SZA = prismaL1_SZA_read(L1_file)

    # Read meanWV from L2C Water Vapor Map product
    meanWV, PRS_L2C_WVM, latitude_WVM, longitude_WVM = prismaL2C_WV_read(L2C_file)

    # Read bounding box from PRISMA L2C file
    bbox = prismaL2C_bbox_read(L2C_file)

    # Analyze the DEM based on the bounding box from PRISMA
    mean_elevation = mean_elev_fromDEM(dem_file, bbox)

    # Extract the file name without the extension for output files
    _, full_filename = os.path.split(L1_file)
    filename_without_extension = os.path.splitext(full_filename)[0]

    # Define output filenames in the specified output directory
    target_spectra_export_name = os.path.join(output_dir, f"{filename_without_extension}_CH4_target_PRISMA_conv.npy")
    mf_output_file = os.path.join(output_dir, f"{filename_without_extension}_MF.tif")
    concentration_output_file = os.path.join(output_dir,  f"{filename_without_extension}_MF_concentration.tif")
    rgb_output_file = os.path.join(output_dir, f"{filename_without_extension}_rgb.tif")
    classified_output_file = os.path.join(output_dir, f"{filename_without_extension}_classified.tif")

    # Data extraction and processing
    rads_array, cw_array, fwhm_array, rgb_image, vnir_cube_bip, swir_cube_bip, latitude_vnir, longitude_vnir, latitude_swir, longitude_swir = prisma_read(L1_file)
    classified_image = k_means_hyperspectral(rads_array, k)
    print("k-means classification completed")
    #plot_classified_image(classified_image)
    rads_array_subselection = rads_array[:, :, 182:214]
    mean_radiance, covariance_matrices = calculate_statistics(rads_array_subselection, classified_image, k)

    # Target spectrum calculation and matched filter application
    concentrations = [0.0, 1000, 2000, 4000, 8000, 16000, 32000, 64000]
    simRads_array, simWave_array = generate_library(concentrations, lut_file, zenith=SZA, sensor=120, ground=mean_elevation, water=meanWV, order=1)
    print("simulated radiance spectrum correctly generated from LUT for 0, 1000, 2000, 4000, 8000, 16000, 32000, 64000 [ppmm] of column CH4 enhancements.")

    # target = generate_template_from_bands(cw_array, fwhm_array, simRads_array, simWave_array, concentrations)
    # np.save(target_spectra_export_name, target)
    #NOTE:The two lines above will be useful for a more efficient implementation that process all band at once instead of relying on the inefficient loop below

    target = None
    for i in range(np.size(cw_array, 0)):
        target_i = generate_template_from_bands(cw_array[i, :], fwhm_array[i, :], simRads_array, simWave_array, concentrations)
        if i == 0:
            target = target_i
        else:
            column_to_add = target_i[:, 1].reshape(-1, 1)
            target = np.concatenate((target, column_to_add), axis=1)

    target_spectra = target[:, 1:]
    # Optional: save to disk target spectra
    np.save(target_spectra_export_name, target_spectra)

    # Selezioniamo solo le bande comprese tra 1500 e 2500 nm
    target_spectra_subselection = target_spectra[182:214, :]

    
    #Matched Filter application
    matched_filter_scores, concentration_map = calculate_matched_filter(rads_array_subselection, classified_image, mean_radiance, covariance_matrices, target_spectra_subselection, k)
    #plot_matched_filter_scores(matched_filter_scores)

    # Save results as GeoTIFF files
    save_as_geotiff_single_band(matched_filter_scores, mf_output_file, latitude_vnir, longitude_vnir)
    save_as_geotiff_single_band(concentration_map , concentration_output_file, latitude_vnir, longitude_vnir)
    save_as_geotiff_rgb(rgb_image, rgb_output_file, latitude_vnir, longitude_vnir)
    save_as_geotiff_single_band(classified_image, classified_output_file, latitude_vnir, longitude_vnir)

    print("output files correctly generated")




##################################################################################################################################################################################
##################################################################################################################################################################################

# This section with argparse must be active if we want the code to run from command line

def parse_args():
    parser = argparse.ArgumentParser(description="Process hyperspectral data.")
    parser.add_argument("L1_file", type=str, help="Input file L1")
    parser.add_argument("L2C_file", type=str, help="Input file L2C")
    parser.add_argument("dem_file", type=str, help="DEM file path")
    parser.add_argument("lut_file", type=str, help="LUT file path")
    parser.add_argument("output_dir", type=str, help="Output directory")
    parser.add_argument("-k", "--clusters", type=int, default=30, help="Number of clusters")
    return parser.parse_args()

def main():
    args = parse_args()
    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Here you would call your main processing function
    ch4_detection(args.L1_file, args.L2C_file, args.dem_file, args.lut_file, args.output_dir, args.clusters)


if __name__ == "__main__":
    main()


##################################################################################################################################################################################
##################################################################################################################################################################################



# # This section with argparse must be active if we want the code to run from editor

# def main(L1_file, L2C_file, dem_file, lut_file, output_dir, k=30):
#     # Ensure output directory exists
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
    
#     # Call the main processing function
#     ch4_detection(L1_file, L2C_file, dem_file, lut_file, output_dir, k)

# #Run in windows
# if __name__ == "__main__":
#     L1_file = r"D:\CLEAR_UP\CH4_detection\img_PRISMA\Turkmenistan\20200807\PRS_L1_STD_OFFL_20200807071742_20200807071747_0001.he5"
#     L2C_file = r"D:\CLEAR_UP\CH4_detection\img_PRISMA\Turkmenistan\20200807\PRS_L2C_STD_20200807071742_20200807071747_0001.he5"
#     dem_file = "D:/CLEAR_UP/CH4_detection/Matched_filter_approach/codici/CTMF/DEM_1Km/srtm30plus_v11_land.nc"
#     lut_file = "D:/CLEAR_UP/CH4_detection/Matched_filter_approach/codici/CTMF/LUTs/dataset_ch4_full.hdf5"
#     output_dir = r"D:\CLEAR_UP\CH4_detection\img_PRISMA\Turkmenistan\20200807\MF_out_trial2"
#     k = 10
#     main(L1_file, L2C_file, dem_file, lut_file, output_dir, k=k)

