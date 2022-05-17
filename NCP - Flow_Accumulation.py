#!/usr/bin/env python
# coding: utf-8

# In[1]:


import logging
import os

import numpy
import pygeoprocessing
import pygeoprocessing.routing
from osgeo import gdal


# In[2]:


import numpy
import pygeoprocessing
from osgeo import gdal

FLOAT32_NODATA = float(numpy.finfo(numpy.float32).min)


# In[3]:


def _convert_to_from_density(source_raster_path, target_raster_path,
                             direction='to_density'):
    """Convert a raster to/from counts/pixel and counts/unit area.
    Args:
        source_raster_path (string): The path to a raster containing units that
            need to be converted.
        target_raster_path (string): The path to where the target raster with
            converted units should be stored.
        direction='to_density' (string): The direction of unit conversion.  If
            'to_density', then the units of ``source_raster_path`` must be in
            counts per pixel and will be converted to counts per square meter.
            If 'from_density', then the units of ``source_raster_path`` must be
            in counts per square meter and will be converted to counts per
            pixel.
    Returns:
        ``None``
    """
    source_raster_info = pygeoprocessing.get_raster_info(source_raster_path)
    source_nodata = source_raster_info['nodata'][0]

    # Calculate the per-pixel area based on the latitude.
    _, miny, _, maxy = source_raster_info['bounding_box']
    pixel_area_in_m2_by_latitude = (
        pygeoprocessing._create_latitude_m2_area_column(
            miny, maxy, source_raster_info['raster_size'][1]))

    def _convert(array, pixel_area):
        out_array = numpy.full(array.shape, FLOAT32_NODATA,
                               dtype=numpy.float32)
        valid_mask = slice(None)
        if source_nodata is not None:
            valid_mask = ~numpy.isclose(array, source_nodata)

        if direction == 'to_density':
            out_array[valid_mask] = array[valid_mask] / pixel_area[valid_mask]
        elif direction == 'from_density':
            out_array[valid_mask] = array[valid_mask] * pixel_area[valid_mask]
        else:
            raise AssertionError(f'Invalid direction: {direction}')
        return out_array

    pygeoprocessing.raster_calculator(
        [(source_raster_path, 1), pixel_area_in_m2_by_latitude],
        _convert, target_raster_path, gdal.GDT_Float32, FLOAT32_NODATA)


# In[5]:


if __name__ == '__main__':
    population_raster_path = 'n_load.tif'
    pop_per_sq_m = 'population_per_sq_m.tif'
    _convert_to_from_density(population_raster_path, pop_per_sq_m,
                             'to_density')

    # call align_and_resize_raster_stack here using pop_per_sq_m as a source
    # path.
    aligned_pop_per_sq_m = 'aligned_pop_per_sq_m.tif'
    pygeoprocessing.align_and_resize_raster_stack(
        base_raster_path_list=[pop_per_sq_m],
        target_raster_path_list=[aligned_pop_per_sq_m],
        resample_method_list=['bilinear'],
        target_pixel_size=pygeoprocessing.get_raster_info(population_raster_path)['pixel_size'],
        bounding_box_mode='intersection'
    )

    # Once pop-per-sq-m has been aligned with the stack, convert it back to
    # pop-per-pixel.
    pop_per_pixel = 'aligned_population_per_pixel.tif'
    _convert_to_from_density(aligned_pop_per_sq_m, pop_per_pixel,
                             'from_density')


# In[ ]:


def convert_taudem_flow_dir_to_pygeoprocessing(
        taudem_flow_dir_raster_path, target_flow_dir_raster_path):
    source_nodata = pygeoprocessing.get_raster_info(
        taudem_flow_dir_raster_path)['nodata'][0]
    dest_nodata = 128  # Matches pygeoprocessing-produced flow dir nodata.

    def _convert(taudem_flow_dir_array):
        dest_array = numpy.full(taudem_flow_dir_array.shape, dest_nodata,
                                dtype=numpy.uint8)
        valid_pixels = (taudem_flow_dir_array != source_nodata)
        dest_array[valid_pixels] = taudem_flow_dir_array[valid_pixels] - 1
        return dest_array

    pygeoprocessing.raster_calculator(
        [(taudem_flow_dir_raster_path, 1)], _convert,
        target_flow_dir_raster_path, gdal.GDT_Byte, dest_nodata)


# In[3]:


if __name__ == '__main__':
    convert_taudem_flow_dir_to_pygeoprocessing(
        'TAUDEM_d8.tif',
        'pygeo_tif.tif')


# In[4]:


logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


# In[6]:


def doit(dem_path, flow_dir_weights,flow_d8_dir, workspace):
    LOGGER.info('Starting script')
    if not os.path.exists(workspace):
        os.makedirs(workspace)

    aligned_dem_path = os.path.join(workspace, 'aligned_dem.tif')
    aligned_weights_path = os.path.join(workspace, 'aligned_weights.tif')
    aligned_flow_d8_dir = os.path.join(workspace, 'aligned_flow_dir_d8.tif' )

    # Align the DEM and weights raster
    pygeoprocessing.align_and_resize_raster_stack(
        [dem_path, flow_dir_weights, flow_d8_dir],
        [aligned_dem_path, aligned_weights_path,aligned_flow_d8_dir ],
        ['bilinear', 'bilinear','bilinear'],
        pygeoprocessing.get_raster_info(dem_path)['pixel_size'],
        bounding_box_mode='intersection')

    # Flow direction
    #flow_direction_path = os.path.join(workspace, 'pygeoprocessing_d8_flow_dir.tif')
    #pygeoprocessing.routing.flow_dir_d8(
     #   (aligned_dem_path, 1), flow_d8_dir)

    # Flow accumulation with weights.
    flow_accumulation_path = os.path.join(workspace, 'flow_accumulation.tif')
    pygeoprocessing.routing.flow_accumulation_d8(
        (aligned_flow_d8_dir, 1), flow_accumulation_path,
        (aligned_weights_path, 1))


# In[7]:


if __name__ == '__main__':
    doit('my_dem_projected.tif', 'aligned_population_per_pixel.tif','pygeo_tif.tif', 'NCP_n_load')


# In[ ]:




