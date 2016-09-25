import properties

location = properties.Vector2(
    "interferogram location (bottom N, left E)",
    required=True
)

location_UTM_zone = properties.Integer(
    "UTM zone",
    required=True
)

interferogram_shape = properties.Array(
    "number of pixels in the interferogram",
    shape=(2,),
    dtype=int,
    required=True
)

interferogram_pixel_size = properties.Array(
    "Size of each pixel (northing, easting)",
    shape=(2,),
    dtype=float,
    required=True
)

interferogram_ref = properties.Vector2(
    "interferogram reference",
    required=True
)

interferogram_ref_incidence = properties.Float(
    "Incidence angle"
)

interferogram_scaling = properties.Float(
    "Scaling of the interferogram",
    default=1.0
)

satellite_name = properties.String("Name of the satelite.")
satellite_fringe_interval = properties.Float(
    "Fringe interval",
    default=0.028333
)

satellite_azimuth = properties.Float("satellite_azimuth")
satellite_altitude = properties.Float("satellite_altitude")

local_rigidity = properties.Float(
    "Local rigidity",
    default=3e10
)

local_earth_radius = properties.Float(
    "Earth radius",
    default=6371000.
)

interferogram_date1 = properties.DateTime(
    "interferogram_date1",
    required=True
)

interferogram_date2 = properties.DateTime(
    "interferogram_date2",
    required=True
)

interferogram_processed_by = properties.String(
    "interferogram_processed_by",
    required=True
)

interferogram_processed_date = properties.DateTime(
    "interferogram_processed_date",
    required=True
)

interferogram_copyright = properties.String(
    "interferogram_copyright",
    required=True
)

interferogram_data_source = properties.String(
    "interferogram_data_source",
    required=True
)

event_date = properties.DateTime("Date of the earthquake")
event_gcmt_id = properties.String("GCMT ID")
event_name = properties.String("Earthquake name")
event_country = properties.String("Earthquake country")

# data_type_p = properties.String("")
