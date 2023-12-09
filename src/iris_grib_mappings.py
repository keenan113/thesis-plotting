
_level_type_mappings = {
    "meanSea": (101, 255),
    "hybrid": (105, 255),
    "atmosphere": (10, 255),
    "atmosphereSingleLayer": (200, 255),
    "surface": (1, 255),
    "isobaricInhPa": (100, 255),
    "heightAboveGround": (103, 255),
    "depthBelowLandLayer": (106, 106),
    "heightAboveSea": (102, 255),
    "isobaricLayer": (100, 100),
    "nominalTop": (8, 255),
    "heightAboveGroundLayer": (103, 103),
    "tropopause": (7, 255),
    "maxWind": (6, 255),
    "isothermZero": (4, 255),
    "pressureFromGroundLayer": (108, 108),
    "sigmaLayer": (104, 104),
    "sigma": (104, 255),
    "theta": (107, 255),
    "potentialVorticity": (109, 255),
}

_field_keyval_mapping = {
    "air_temperature_isobaric":{
        "discipline":0,
        "category":0,
        "number":0,
        "type_of_level":"isobaricInhPa",
        "product_definition_template":0,
    },
    "specific_humidity_isobaric":{
        "discipline":0,
        "category":1,
        "number":0,
        "type_of_level":"isobaricInhPa",
        "product_definition_template":0,
    },
    "geopotential_height_isobaric":{
        "discipline":0,
        "category":3,
        "number":5,
        "type_of_level":"isobaricInhPa",
        "product_definition_template":0,
    },
    "u_wind_isobaric":{
        "discipline":0,
        "category":2,
        "number":2,
        "type_of_level":"isobaricInhPa",
        "product_definition_template":0,
    },
    "v_wind_isobaric":{
        "discipline":0,
        "category":2,
        "number":3,
        "type_of_level":"isobaricInhPa",
        "product_definition_template":0,
    },
    "latitude":{
        "discipline":0,
        "category":191,
        "number":192,
        "type_of_level":"surface",
        "product_definition_template":0
    },
    "longitude":{
        "discipline":0,
        "category":191,
        "number":193,
        "type_of_level":"surface",
        "product_definition_template":0
    },
    "u_wind_surface":{  #U-component of winds at 10m AGL (m/s)
        "discipline":0,
        "category":2,
        "number":2,
        "type_of_level":"heightAboveGround",
        "product_definition_template":0,
        "level":10
    },
    "v_wind_surface":{  #V-component of winds at 10m AGL (m/s)
        "discipline":0,
        "category":2,
        "number":3,
        "type_of_level":"heightAboveGround",
        "product_definition_template":0,
        "level":10
    },
    "pressure_surface":{    #surface pressure (Pa)
        "discipline":0,
        "category":3,
        "number":0,
        "type_of_level":"surface",
        "product_definition_template":0
    },
    "air_temperature_surface":{ #air temperature 2-meters above ground
        "discipline":0,
        "category":0,
        "number":0,
        "type_of_level":"heightAboveGround",
        "product_definition_template":0,
        "level":2
    },
    "dewpoint_temperature_surface":{ #dew point temperature 2-meters above ground
        "discipline":0,
        "category":0,
        "number":6,
        "type_of_level":"heightAboveGround",
        "product_definition_template":0,
        "level":2
    },
    "terrain_height_surface":{  #Height (above sea-level) of the model terrain
        "discipline":0,
        "category":3,
        "number":5,
        "type_of_level":"surface",
        "product_definition_template":0
    },
    "mslp_surface":{  #Mean Sea Level Pressure
        "discipline":0,
        "category":3,
        "number":1,
        "type_of_level":"meanSea",
        "product_definition_template":0
    },
    "mslp_MAPS_surface":{
        "discipline":0,
        "category":3,
        "number":198,
        "type_of_level":"meanSea",
        "product_definition_template":0
    },
    "composite_reflectivity_surface":{  #Max/Composite Reflectivity
        "discipline":0,
        "category":16,
        "number":196,
        "type_of_level":"atmosphereSingleLayer",
        "product_definition_template":0
    },
#    "total_precipitation_surface":{  #run-total precipitation 
#        "discipline":0,
#        "category":1,
#        "number":8,
#        "type_of_level":"surface",
#        "product_definition_template":8
#    },
}
