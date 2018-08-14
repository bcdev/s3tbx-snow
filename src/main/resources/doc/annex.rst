.. _annex:

=====
Annex
=====

Example of IdePix NetCDF4 product header:
::

    netcdf L2_of_MER_RR__1PNACR20080621_055731_000001512069_00363_32982_0000 {
    dimensions:
        y = 865 ;
        x = 1121 ;
        tp_y = 55 ;
        tp_x = 71 ;
    variables:
        short cloud_classif_flags(y, x) ;
            cloud_classif_flags:coordinates = "lat lon" ;
            cloud_classif_flags:flag_meanings = "F_INVALID F_CLOUD
            F_CLOUD_AMBIGUOUS F_CLOUD_SURE F_CLOUD_BUFFER F_CLOUD_SHADOW
            F_SNOW_ICE F_GLINTRISK F_COASTLINE F_LAND" ;
            cloud_classif_flags:flag_masks = 1s, 2s, 4s, 8s, 16s, 32s, 64s,
            128s, 256s, 512s ;
            cloud_classif_flags:flag_coding_name = "cloud_classif_flags" ;
            cloud_classif_flags:flag_descriptions = "Invalid pixels\tPixels
            which are either cloud_sure or cloud_ambiguous\tSemi transparent
            clouds, or clouds where the detection level is uncertain\tFully
            opaque clouds with full confidence of their detection\tA buffer
            of n pixels around a cloud. n is a user supplied parameter. Applied
            to pixels masked as \'cloud\'\tPixels is affect by a cloud
            shadow\tSnow/ice pixels\tPixels with glint risk\tPixels at a
            coastline\tLand pixels" ;
            cloud_classif_flags:long_name = "" ;
        short radiance_10(y, x) ;
            radiance_10:long_name = "TOA radiance band 10" ;
            radiance_10:units = "mW/(m^2*sr*nm)" ;
            radiance_10:_Unsigned = "true" ;
            radiance_10:scale_factor = 0.00866463407874107 ;
            radiance_10:coordinates = "lat lon" ;
            radiance_10:bandwidth = 7.495f ;
            radiance_10:wavelength = 753.371f ;
            radiance_10:valid_pixel_expression = "!l1_flags.INVALID" ;
            radiance_10:solar_flux = 1227.051f ;
            radiance_10:spectral_band_index = 9.f ;
        short radiance_11(y, x) ;
            radiance_11:long_name = "TOA radiance band 11" ;
            radiance_11:units = "mW/(m^2*sr*nm)" ;
            radiance_11:_Unsigned = "true" ;
            radiance_11:scale_factor = 0.00887294951826334 ;
            radiance_11:coordinates = "lat lon" ;
            radiance_11:bandwidth = 3.744f ;
            radiance_11:wavelength = 761.5081f ;
            radiance_11:valid_pixel_expression = "!l1_flags.INVALID" ;
            radiance_11:solar_flux = 1215.942f ;
            radiance_11:spectral_band_index = 10.f ;
        short detector_index(y, x) ;
            detector_index:coordinates = "lat lon" ;
            detector_index:long_name = "Detector index" ;
        byte l1_flags(y, x) ;
            l1_flags:_Unsigned = "true" ;
            l1_flags:coordinates = "lat lon" ;
            l1_flags:flag_meanings = "COSMETIC DUPLICATED GLINT_RISK SUSPECT
            LAND_OCEAN BRIGHT COASTLINE INVALID" ;
            l1_flags:flag_masks = 1b, 2b, 4b, 8b, 16b, 32b, 64b, -128b ;
            l1_flags:flag_coding_name = "l1_flags" ;
            l1_flags:flag_descriptions = "Pixel is cosmetic\tPixel has been
            duplicated (filled in)\tPixel has glint risk\tPixel is
            suspect\tPixel is over land, not ocean\tPixel is bright\tPixel
            is part of acoastline\tPixel is invalid" ;
            l1_flags:long_name = "Level 1b classification and quality flags" ;
        float latitude(tp_y, tp_x) ;
            latitude:offset_y = 0.5f ;
            latitude:subsampling_x = 16.f ;
            latitude:subsampling_y = 16.f ;
            latitude:offset_x = 0.5f ;
        float longitude(tp_y, tp_x) ;
            longitude:offset_y = 0.5f ;
            longitude:subsampling_x = 16.f ;
            longitude:subsampling_y = 16.f ;
            longitude:offset_x = 0.5f ;
        float dem_alt(tp_y, tp_x) ;
            dem_alt:offset_y = 0.5f ;
            dem_alt:subsampling_x = 16.f ;
            dem_alt:subsampling_y = 16.f ;
            dem_alt:offset_x = 0.5f ;
        float dem_rough(tp_y, tp_x) ;
            dem_rough:offset_y = 0.5f ;
            dem_rough:subsampling_x = 16.f ;
            dem_rough:subsampling_y = 16.f ;
            dem_rough:offset_x = 0.5f ;
        float lat_corr(tp_y, tp_x) ;
            lat_corr:offset_y = 0.5f ;
            lat_corr:subsampling_x = 16.f ;
            lat_corr:subsampling_y = 16.f ;
            lat_corr:offset_x = 0.5f ;
        float lon_corr(tp_y, tp_x) ;
            lon_corr:offset_y = 0.5f ;
            lon_corr:subsampling_x = 16.f ;
            lon_corr:subsampling_y = 16.f ;
            lon_corr:offset_x = 0.5f ;
        float sun_zenith(tp_y, tp_x) ;
            sun_zenith:offset_y = 0.5f ;
            sun_zenith:subsampling_x = 16.f ;
            sun_zenith:subsampling_y = 16.f ;
            sun_zenith:offset_x = 0.5f ;
        float sun_azimuth(tp_y, tp_x) ;
            sun_azimuth:offset_y = 0.5f ;
            sun_azimuth:subsampling_x = 16.f ;
            sun_azimuth:subsampling_y = 16.f ;
            sun_azimuth:offset_x = 0.5f ;
        float view_zenith(tp_y, tp_x) ;
            view_zenith:offset_y = 0.5f ;
            view_zenith:subsampling_x = 16.f ;
            view_zenith:subsampling_y = 16.f ;
            view_zenith:offset_x = 0.5f ;
        float view_azimuth(tp_y, tp_x) ;
            view_azimuth:offset_y = 0.5f ;
            view_azimuth:subsampling_x = 16.f ;
            view_azimuth:subsampling_y = 16.f ;
            view_azimuth:offset_x = 0.5f ;
        float zonal_wind(tp_y, tp_x) ;
            zonal_wind:offset_y = 0.5f ;
            zonal_wind:subsampling_x = 16.f ;
            zonal_wind:subsampling_y = 16.f ;
            zonal_wind:offset_x = 0.5f ;
        float merid_wind(tp_y, tp_x) ;
            merid_wind:offset_y = 0.5f ;
            merid_wind:subsampling_x = 16.f ;
            merid_wind:subsampling_y = 16.f ;
            merid_wind:offset_x = 0.5f ;
        float atm_press(tp_y, tp_x) ;
            atm_press:offset_y = 0.5f ;
            atm_press:subsampling_x = 16.f ;
            atm_press:subsampling_y = 16.f ;
            atm_press:offset_x = 0.5f ;
        float ozone(tp_y, tp_x) ;
            ozone:offset_y = 0.5f ;
            ozone:subsampling_x = 16.f ;
            ozone:subsampling_y = 16.f ;
            ozone:offset_x = 0.5f ;
        float rel_hum(tp_y, tp_x) ;
            rel_hum:offset_y = 0.5f ;
            rel_hum:subsampling_x = 16.f ;
            rel_hum:subsampling_y = 16.f ;
            rel_hum:offset_x = 0.5f ;
        float lat(y, x) ;
            lat:long_name = "latitude coordinate" ;
            lat:standard_name = "latitude" ;
            lat:units = "degrees_north" ;
        float lon(y, x) ;
            lon:long_name = "longitude coordinate" ;
            lon:standard_name = "longitude" ;
            lon:units = "degrees_east" ;
        byte cawa_invalid_mask ;
            cawa_invalid_mask:expression = "cloud_classif_flags.F_INVALID" ;
            cawa_invalid_mask:color = 178, 0, 0, 255 ;
            cawa_invalid_mask:transparency = 0.5 ;
            cawa_invalid_mask:title = "Invalid pixels" ;
        byte cawa_cloud_mask ;
            cawa_cloud_mask:expression = "cloud_classif_flags.F_CLOUD" ;
            cawa_cloud_mask:color = 255, 0, 255, 255 ;
            cawa_cloud_mask:transparency = 0.5 ;
            cawa_cloud_mask:title = "Pixels which are either cloud_sure or
            cloud_ambiguous" ;
        byte cawa_cloud_ambiguous_mask ;
            cawa_cloud_ambiguous_mask:expression =
            "cloud_classif_flags.F_CLOUD_AMBIGUOUS" ;
            cawa_cloud_ambiguous_mask:color = 255, 255, 0, 255 ;
            cawa_cloud_ambiguous_mask:transparency = 0.5 ;
            cawa_cloud_ambiguous_mask:title = "Semi transparent clouds, or
            clouds where the detection level  is uncertain" ;
        byte cawa_cloud_sure_mask ;
            cawa_cloud_sure_mask:expression =
            "cloud_classif_flags.F_CLOUD_SURE" ;
            cawa_cloud_sure_mask:color = 255, 0, 0, 255 ;
            cawa_cloud_sure_mask:transparency = 0.5 ;
            cawa_cloud_sure_mask:title = "Fully opaque clouds with full
            confidence of their detection" ;
        byte cawa_cloud_buffer_mask ;
            cawa_cloud_buffer_mask:expression =
            "cloud_classif_flags.F_CLOUD_BUFFER" ;
            cawa_cloud_buffer_mask:color = 255, 200, 0, 255 ;
            cawa_cloud_buffer_mask:transparency = 0.5 ;
            cawa_cloud_buffer_mask:title = "A buffer of n pixels around a
            cloud. n is a user supplied parameter.
            Applied to pixels masked as \'cloud\'" ;
        byte cawa_cloud_shadow_mask ;
            cawa_cloud_shadow_mask:expression =
            "cloud_classif_flags.F_CLOUD_SHADOW" ;
            cawa_cloud_shadow_mask:color = 178, 0, 0, 255 ;
            cawa_cloud_shadow_mask:transparency = 0.5 ;
            cawa_cloud_shadow_mask:title =
            "Pixels is affect by a cloud shadow" ;
        byte cawa_snow_ice_mask ;
            cawa_snow_ice_mask:expression = "cloud_classif_flags.F_SNOW_ICE" ;
            cawa_snow_ice_mask:color = 0, 255, 255, 255 ;
            cawa_snow_ice_mask:transparency = 0.5 ;
            cawa_snow_ice_mask:title = "Snow/ice pixels" ;
        byte cawa_glint_risk_mask ;
            cawa_glint_risk_mask:expression =
            "cloud_classif_flags.F_GLINTRISK" ;
            cawa_glint_risk_mask:color = 255, 175, 175, 255 ;
            cawa_glint_risk_mask:transparency = 0.5 ;
            cawa_glint_risk_mask:title = "Pixels with glint risk" ;
        byte cawa_coastline_mask ;
            cawa_coastline_mask:expression = "cloud_classif_flags.F_COASTLINE" ;
            cawa_coastline_mask:color = 0, 178, 0, 255 ;
            cawa_coastline_mask:transparency = 0.5 ;
            cawa_coastline_mask:title = "Pixels at a coastline" ;
        byte cawa_land_mask ;
            cawa_land_mask:expression = "cloud_classif_flags.F_LAND" ;
            cawa_land_mask:color = 0, 255, 0, 255 ;
            cawa_land_mask:transparency = 0.5 ;
            cawa_land_mask:title = "Land pixels" ;
        byte coastline_mask ;
            coastline_mask:expression = "l1_flags.COASTLINE" ;
            coastline_mask:color = 0, 255, 0, 255 ;
            coastline_mask:transparency = 0. ;
            coastline_mask:title = "Pixel is part of a coastline" ;
        byte land_mask ;
            land_mask:expression = "l1_flags.LAND_OCEAN" ;
            land_mask:color = 51, 153, 0, 255 ;
            land_mask:transparency = 0.75 ;
            land_mask:title = "Pixel is over land, not ocean" ;
        byte water_mask ;
            water_mask:expression = "NOT l1_flags.LAND_OCEAN" ;
            water_mask:color = 153, 153, 255, 255 ;
            water_mask:transparency = 0.75 ;
            water_mask:title = "Not Pixel is over land, not ocean" ;
        byte cosmetic_mask ;
            cosmetic_mask:expression = "l1_flags.COSMETIC" ;
            cosmetic_mask:color = 204, 153, 255, 255 ;
            cosmetic_mask:transparency = 0.5 ;
            cosmetic_mask:title = "Pixel is cosmetic" ;
        byte duplicated_mask ;
            duplicated_mask:expression = "l1_flags.DUPLICATED" ;
            duplicated_mask:color = 255, 200, 0, 255 ;
            duplicated_mask:transparency = 0.5 ;
            duplicated_mask:title = "Pixel has been duplicated (filled in)" ;
        byte glint_risk_mask ;
            glint_risk_mask:expression = "l1_flags.GLINT_RISK" ;
            glint_risk_mask:color = 255, 0, 255, 255 ;
            glint_risk_mask:transparency = 0.5 ;
            glint_risk_mask:title = "Pixel has glint risk" ;
        byte suspect_mask ;
            suspect_mask:expression = "l1_flags.SUSPECT" ;
            suspect_mask:color = 204, 102, 255, 255 ;
            suspect_mask:transparency = 0.5 ;
            suspect_mask:title = "Pixel is suspect" ;
        byte bright_mask ;
            bright_mask:expression = "l1_flags.BRIGHT" ;
            bright_mask:color = 255, 255, 0, 255 ;
            bright_mask:transparency = 0.5 ;
            bright_mask:title = "Pixel is bright" ;
        byte invalid_mask ;
            invalid_mask:expression = "l1_flags.INVALID" ;
            invalid_mask:color = 255, 0, 0, 255 ;
            invalid_mask:transparency = 0. ;
            invalid_mask:title = "Pixel is invalid" ;

    // global attributes:
            :Conventions = "CF-1.4" ;
            :TileSize = "16:1121" ;
            :product_type = "mergedClassif" ;
            :metadata_profile = "beam" ;
            :metadata_version = "0.5" ;
            :auto_grouping = "radiance:rho_toa" ;
            :tiepoint_coordinates = "longitude latitude" ;
            :start_date = "21-JUN-2008 05:57:31.155941" ;
            :stop_date = "21-JUN-2008 06:00:03.209572" ;
    }


Example of CAWA TCWV product header:
::

    netcdf L2_of_L2_of_MER_RR__1PNUPA20060102_141100_000026182043_00497_20090_7596 {
    dimensions:
        y = 14881 ;
        x = 1121 ;
        tp_y = 931 ;
        tp_x = 71 ;
    variables:
        float tcwv(y, x) ;
            tcwv:units = "mm" ;
            tcwv:_FillValue = -999.f ;
            tcwv:long_name = "Total column of water vapour" ;
        byte tcwv_flags(y, x) ;
            tcwv_flags:units = "1" ;
            tcwv_flags:long_name = "TCWV flags band" ;
        short cloud_classif_flags(y, x) ;
            cloud_classif_flags:units = "1" ;
            cloud_classif_flags:flag_meanings = "F_INVALID F_CLOUD
            F_CLOUD_AMBIGUOUS F_CLOUD_SURE F_CLOUD_BUFFER F_CLOUD_SHADOW
            F_SNOW_ICE F_GLINTRISK F_COASTLINE F_LAND" ;
            cloud_classif_flags:flag_masks = 1s, 2s, 4s, 8s, 16s, 32s, 64s,
            128s, 256s, 512s ;
            cloud_classif_flags:flag_coding_name = "cloud_classif_flags" ;
            cloud_classif_flags:flag_descriptions = "Invalid pixels\tPixels
            which are either cloud_sure or cloud_ambiguous\tSemi transparent
            clouds, or clouds where the detection level is uncertain\tFully
            opaque clouds with full confidence of their detection\tA buffer
            of n pixels around a cloud. n is a user supplied parameter. Applied
            to pixels masked as \'cloud\'\tPixels is affect by a cloud
            shadow\tSnow/ice pixels\tPixels with glint risk\tPixels at a
            coastline\tLand pixels" ;
            cloud_classif_flags:long_name = "" ;
        float latitude(tp_y, tp_x) ;
            latitude:offset_y = 0.5 ;
            latitude:subsampling_x = 16. ;
            latitude:subsampling_y = 16. ;
            latitude:units = "degree" ;
            latitude:standard_name = "latitude" ;
            latitude:offset_x = 0.5 ;
        float longitude(tp_y, tp_x) ;
            longitude:offset_y = 0.5 ;
            longitude:subsampling_x = 16. ;
            longitude:subsampling_y = 16. ;
            longitude:units = "degree" ;
            longitude:standard_name = "longitude" ;
            longitude:offset_x = 0.5 ;
        byte cawa_invalid_mask ;
            cawa_invalid_mask:description = "Invalid pixels" ;
            cawa_invalid_mask:expression = "cloud_classif_flags.F_INVALID" ;
            cawa_invalid_mask:color = 178, 0, 0, 255 ;
            cawa_invalid_mask:transparency = 0.5 ;
            cawa_invalid_mask:long_name = "cawa_invalid" ;
        byte cawa_cloud_mask ;
            cawa_cloud_mask:description = "Pixels which are either cloud_sure
            or cloud_ambiguous" ;
            cawa_cloud_mask:expression = "cloud_classif_flags.F_CLOUD" ;
            cawa_cloud_mask:color = 255, 0, 255, 255 ;
            cawa_cloud_mask:transparency = 0.5 ;
            cawa_cloud_mask:long_name = "cawa_cloud" ;
        byte cawa_cloud_ambiguous_mask ;
            cawa_cloud_ambiguous_mask:description = "Semi transparent clouds,
            or clouds where the detection level is uncertain" ;
            cawa_cloud_ambiguous_mask:expression =
            "cloud_classif_flags.F_CLOUD_AMBIGUOUS" ;
            cawa_cloud_ambiguous_mask:color = 255, 255, 0, 255 ;
            cawa_cloud_ambiguous_mask:transparency = 0.5 ;
            cawa_cloud_ambiguous_mask:long_name = "cawa_cloud_ambiguous" ;
        byte cawa_cloud_sure_mask ;
            cawa_cloud_sure_mask:description =
            "Fully opaque clouds with full confidence of their detection" ;
            cawa_cloud_sure_mask:expression = "cloud_classif_flags.F_CLOUD_SURE" ;
            cawa_cloud_sure_mask:color = 255, 0, 0, 255 ;
            cawa_cloud_sure_mask:transparency = 0.5 ;
            cawa_cloud_sure_mask:long_name = "cawa_cloud_sure" ;
        byte cawa_cloud_buffer_mask ;
            cawa_cloud_buffer_mask:description = "A buffer of n pixels around
            a cloud. n is a user supplied parameter. Applied to pixels masked
            as \'cloud\'" ;
            cawa_cloud_buffer_mask:expression =
            "cloud_classif_flags.F_CLOUD_BUFFER" ;
            cawa_cloud_buffer_mask:color = 255, 200, 0, 255 ;
            cawa_cloud_buffer_mask:transparency = 0.5 ;
            cawa_cloud_buffer_mask:long_name = "cawa_cloud_buffer" ;
        byte cawa_cloud_shadow_mask ;
            cawa_cloud_shadow_mask:description = "Pixels is affect by a
            cloud shadow" ;
            cawa_cloud_shadow_mask:expression =
            "cloud_classif_flags.F_CLOUD_SHADOW" ;
            cawa_cloud_shadow_mask:color = 178, 0, 0, 255 ;
            cawa_cloud_shadow_mask:transparency = 0.5 ;
            cawa_cloud_shadow_mask:long_name = "cawa_cloud_shadow" ;
        byte cawa_snow_ice_mask ;
            cawa_snow_ice_mask:description = "Snow/ice pixels" ;
            cawa_snow_ice_mask:expression = "cloud_classif_flags.F_SNOW_ICE" ;
            cawa_snow_ice_mask:color = 0, 255, 255, 255 ;
            cawa_snow_ice_mask:transparency = 0.5 ;
            cawa_snow_ice_mask:long_name = "cawa_snow_ice" ;
        byte cawa_glint_risk_mask ;
            cawa_glint_risk_mask:description = "Pixels with glint risk" ;
            cawa_glint_risk_mask:expression = "cloud_classif_flags.F_GLINTRISK" ;
            cawa_glint_risk_mask:color = 255, 175, 175, 255 ;
            cawa_glint_risk_mask:transparency = 0.5 ;
            cawa_glint_risk_mask:long_name = "cawa_glint_risk" ;
        byte cawa_coastline_mask ;
            cawa_coastline_mask:description = "Pixels at a coastline" ;
            cawa_coastline_mask:expression = "cloud_classif_flags.F_COASTLINE" ;
            cawa_coastline_mask:color = 0, 178, 0, 255 ;
            cawa_coastline_mask:transparency = 0.5 ;
            cawa_coastline_mask:long_name = "cawa_coastline" ;
        byte cawa_land_mask ;
            cawa_land_mask:description = "Land pixels" ;
            cawa_land_mask:expression = "cloud_classif_flags.F_LAND" ;
            cawa_land_mask:color = 0, 255, 0, 255 ;
            cawa_land_mask:transparency = 0.5 ;
            cawa_land_mask:long_name = "cawa_land" ;

    // global attributes:
            :Conventions = "CF-1.4" ;
            :title = "CAWA TCWV product" ;
            :product_type = "CAWA TCWV" ;
            :start_date = "02-JAN-2006 14:11:00.727666" ;
            :stop_date = "02-JAN-2006 14:54:39.429106" ;
            :TileSize = "64:1121" ;
            :metadata_profile = "beam" ;
            :metadata_version = "0.5" ;
            :tiepoint_coordinates = "longitude latitude" ;
    }



Example of CAWA CTP product header:
::

    netcdf L2_of_L2_of_MER_RR__1PNUPA20050701_072830_000026412038_00350_17438_5743 {
    dimensions:
        y = 15009 ;
        x = 1121 ;
        tp_y = 939 ;
        tp_x = 71 ;
    variables:
        float ctp(y, x) ;
            ctp:units = "hPa" ;
            ctp:_FillValue = -999.f ;
            ctp:long_name = "Cloud Top Pressure" ;
        byte ctp_flags(y, x) ;
            ctp_flags:units = "1" ;
            ctp_flags:long_name = "CTP flags band" ;
        short cloud_classif_flags(y, x) ;
            cloud_classif_flags:units = "1" ;
            cloud_classif_flags:flag_meanings = "F_INVALID F_CLOUD
            F_CLOUD_AMBIGUOUS F_CLOUD_SURE F_CLOUD_BUFFER F_CLOUD_SHADOW
            F_SNOW_ICE F_GLINTRISK F_COASTLINE F_LAND" ;
            cloud_classif_flags:flag_masks = 1s, 2s, 4s, 8s, 16s, 32s, 64s,
            128s, 256s, 512s ;
            cloud_classif_flags:flag_coding_name = "cloud_classif_flags" ;
            cloud_classif_flags:flag_descriptions = "Invalid pixels\tPixels
            which are either cloud_sure or cloud_ambiguous\tSemi transparent
            clouds, or clouds where the detection level is uncertain\tFully
            opaque clouds with full confidence of their detection\tA buffer
            of n pixels around a cloud. n is a user supplied parameter.
            Applied to pixels masked as \'cloud\'\tPixels is affect by a cloud
            shadow\tSnow/ice pixels\tPixels with glint risk\tPixels at a
            coastline\tLand pixels" ;
            cloud_classif_flags:long_name = "" ;
        float latitude(tp_y, tp_x) ;
            latitude:offset_y = 0.5 ;
            latitude:subsampling_x = 16. ;
            latitude:subsampling_y = 16. ;
            latitude:units = "degree" ;
            latitude:standard_name = "latitude" ;
            latitude:offset_x = 0.5 ;
        float longitude(tp_y, tp_x) ;
            longitude:offset_y = 0.5 ;
            longitude:subsampling_x = 16. ;
            longitude:subsampling_y = 16. ;
            longitude:units = "degree" ;
            longitude:standard_name = "longitude" ;
            longitude:offset_x = 0.5 ;
        byte cawa_invalid_mask ;
            cawa_invalid_mask:description = "Invalid pixels" ;
            cawa_invalid_mask:expression = "cloud_classif_flags.F_INVALID" ;
            cawa_invalid_mask:color = 178, 0, 0, 255 ;
            cawa_invalid_mask:transparency = 0.5 ;
            cawa_invalid_mask:long_name = "cawa_invalid" ;
        byte cawa_cloud_mask ;
            cawa_cloud_mask:description = "Pixels which are either cloud_sure
            or cloud_ambiguous" ;
            cawa_cloud_mask:expression = "cloud_classif_flags.F_CLOUD" ;
            cawa_cloud_mask:color = 255, 0, 255, 255 ;
            cawa_cloud_mask:transparency = 0.5 ;
            cawa_cloud_mask:long_name = "cawa_cloud" ;
        byte cawa_cloud_ambiguous_mask ;
            cawa_cloud_ambiguous_mask:description = "Semi transparent clouds,
            or clouds where the detection level is uncertain" ;
            cawa_cloud_ambiguous_mask:expression =
            "cloud_classif_flags.F_CLOUD_AMBIGUOUS" ;
            cawa_cloud_ambiguous_mask:color = 255, 255, 0, 255 ;
            cawa_cloud_ambiguous_mask:transparency = 0.5 ;
            cawa_cloud_ambiguous_mask:long_name = "cawa_cloud_ambiguous" ;
        byte cawa_cloud_sure_mask ;
            cawa_cloud_sure_mask:description = "Fully opaque clouds with full
            confidence of their detection" ;
            cawa_cloud_sure_mask:expression = "
            cloud_classif_flags.F_CLOUD_SURE" ;
            cawa_cloud_sure_mask:color = 255, 0, 0, 255 ;
            cawa_cloud_sure_mask:transparency = 0.5 ;
            cawa_cloud_sure_mask:long_name = "cawa_cloud_sure" ;
        byte cawa_cloud_buffer_mask ;
            cawa_cloud_buffer_mask:description = "A buffer of n pixels around
            a cloud. n is a user supplied parameter. Applied to pixels masked
            as \'cloud\'" ;
            cawa_cloud_buffer_mask:expression =
            "cloud_classif_flags.F_CLOUD_BUFFER" ;
            cawa_cloud_buffer_mask:color = 255, 200, 0, 255 ;
            cawa_cloud_buffer_mask:transparency = 0.5 ;
            cawa_cloud_buffer_mask:long_name = "cawa_cloud_buffer" ;
        byte cawa_cloud_shadow_mask ;
            cawa_cloud_shadow_mask:description = "Pixels is affect by a
            cloud shadow" ;
            cawa_cloud_shadow_mask:expression =
            "cloud_classif_flags.F_CLOUD_SHADOW" ;
            cawa_cloud_shadow_mask:color = 178, 0, 0, 255 ;
            cawa_cloud_shadow_mask:transparency = 0.5 ;
            cawa_cloud_shadow_mask:long_name = "cawa_cloud_shadow" ;
        byte cawa_snow_ice_mask ;
            cawa_snow_ice_mask:description = "Snow/ice pixels" ;
            cawa_snow_ice_mask:expression = "cloud_classif_flags.F_SNOW_ICE" ;
            cawa_snow_ice_mask:color = 0, 255, 255, 255 ;
            cawa_snow_ice_mask:transparency = 0.5 ;
            cawa_snow_ice_mask:long_name = "cawa_snow_ice" ;
        byte cawa_glint_risk_mask ;
            cawa_glint_risk_mask:description = "Pixels with glint risk" ;
            cawa_glint_risk_mask:expression =
            "cloud_classif_flags.F_GLINTRISK" ;
            cawa_glint_risk_mask:color = 255, 175, 175, 255 ;
            cawa_glint_risk_mask:transparency = 0.5 ;
            cawa_glint_risk_mask:long_name = "cawa_glint_risk" ;
        byte cawa_coastline_mask ;
            cawa_coastline_mask:description = "Pixels at a coastline" ;
            cawa_coastline_mask:expression =
            "cloud_classif_flags.F_COASTLINE" ;
            cawa_coastline_mask:color = 0, 178, 0, 255 ;
            cawa_coastline_mask:transparency = 0.5 ;
            cawa_coastline_mask:long_name = "cawa_coastline" ;
        byte cawa_land_mask ;
            cawa_land_mask:description = "Land pixels" ;
            cawa_land_mask:expression = "cloud_classif_flags.F_LAND" ;
            cawa_land_mask:color = 0, 255, 0, 255 ;
            cawa_land_mask:transparency = 0.5 ;
            cawa_land_mask:long_name = "cawa_land" ;

    // global attributes:
            :Conventions = "CF-1.4" ;
            :title = "CAWA product" ;
            :product_type = "CAWA CTP" ;
            :start_date = "01-JUL-2005 07:28:30.062937" ;
            :stop_date = "01-JUL-2005 08:12:31.290841" ;
            :TileSize = "64:1121" ;
            :metadata_profile = "beam" ;
            :metadata_version = "0.5" ;
            :tiepoint_coordinates = "longitude latitude" ;
    }
