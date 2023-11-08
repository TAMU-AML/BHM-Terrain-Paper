## Dataset

The datasets used in this paper comprises of individual turbine operations data as well as terrain information data. The code will automatically read the appropriate dataset for each bash script. The code files assumes that the data is in `data/` folder and have the same name as given below:

| Files | Description |
|-------|-------------|
| `Turbine1_2017.csv` to `Turbine66_2017.csv` | These files contain 10-min frequency wind power curve data for each of the 66 turbines used in the paper. |
|||
| `thinned_power.csv`, `thinned_wind_speed.csv`, `thinned_temperature.csv` | These files contains thinned data using the thinning procedure described in the paper. Each column corresponds to a single turbine's data in each file. |
|||
| `weightedTerrainData.csv` | This files contains the terrain variables values for each of the turbine weighted by the frequency of data falling in a wind direction sector as described in the paper.|
|||
| `location.csv`| This file contains the location (latitude and longitude) of each turbine. The coordinates have been shifted by an unknown constant to preserve anonymity of the location. |
|||
| `terrain_level_agg.csv` | This file contains the categorical level of the terrain for each turbine as described in the paper. |

