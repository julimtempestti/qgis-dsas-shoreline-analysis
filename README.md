# qgis-dsas-shoreline-analysis
# DSAS-style Shoreline Change Analysis for QGIS

A Python script for QGIS that implements the core shoreline change analysis methodologies of the USGS Digital Shoreline Analysis System (DSAS v5.1), allowing coastal change assessment directly within QGIS without installing DSAS.

## Overview

This tool generates transects from a user-defined baseline, intersects them with historical shoreline positions, and computes shoreline change statistics following the methodologies described in:

> Himmelstoss, E.A., Henderson, R.E., Kratzmann, M.G., and Farris, A.S., 2021. *Digital Shoreline Analysis System (DSAS) version 5.1 User Guide*. U.S. Geological Survey Open-File Report 2021–1091.

The script was developed as a standalone PyQGIS implementation to facilitate shoreline change analysis in environments where DSAS is not available.

## Features

### Automatic transect generation

* Configurable transect spacing
* Configurable transect length
* Smoothed perpendicular orientation
* Automatic land/sea side detection
* Handling of multiple shoreline intersections

### Shoreline change statistics

For each transect, the script calculates:

| Metric | Description                         |
| ------ | ----------------------------------- |
| EPR    | End Point Rate (m/year)             |
| LRR    | Linear Regression Rate (m/year)     |
| WLR    | Weighted Linear Regression (m/year) |
| NSM    | Net Shoreline Movement (m)          |
| SCE    | Shoreline Change Envelope (m)       |

### Forecasting

* Shoreline position projection using LRR
* 95% prediction intervals
* Spatial smoothing of forecast boundaries
* Forecast uncertainty polygon generation

### Visualization

* Automatic DSAS-style symbology
* Red-to-blue erosion/accretion color ramp
* Percentile-based classification
* Ready-to-use layers in QGIS

---

## Inputs

### Baseline

A polyline located landward of the shoreline and approximately parallel to the coast.

### Shorelines

A vector layer containing all historical shoreline positions.

Each feature must contain a date field using one of the following formats:

* DD/MM/YYYY
* YYYY-MM-DD
* YYYY/MM/DD
* DD-MM-YYYY

### Positional uncertainty

A default positional uncertainty can be assigned to all shorelines, with optional date-specific uncertainties for individual datasets or sensors.

---

## Outputs

The script automatically generates:

| Layer              | Description                     |
| ------------------ | ------------------------------- |
| DSAS_EPR           | End Point Rate                  |
| DSAS_LRR           | Linear Regression Rate          |
| DSAS_WLR           | Weighted Linear Regression      |
| DSAS_NSM           | Net Shoreline Movement          |
| DSAS_SCE           | Shoreline Change Envelope       |
| DSAS_Forecast_Lin  | Forecast shoreline position     |
| DSAS_Forecast_Zona | Forecast uncertainty zone (95%) |

Outputs can remain as temporary QGIS layers or be exported directly as shapefiles.

---

## Requirements

* QGIS 3.x
* Python environment included with QGIS
* Shorelines and baseline in a projected CRS with metric units

No additional plugins or external dependencies are required.

---

## Installation

1. Download the script.
2. Open QGIS.
3. Load your baseline and shoreline layers.
4. Open:

```
Plugins → Python Console → Show Editor
```

5. Open the script.
6. Configure the parameters in the **CONFIGURATION** section.
7. Run the script.

---

## Configuration Parameters

Main user-configurable settings include:

* Baseline layer name
* Shoreline layer name
* Date field
* Transect spacing
* Transect length
* Smoothing distance
* Land/sea orientation
* Positional uncertainty
* Forecast horizon
* Output directory

No code modification is required outside the configuration section.

---

## Typical Applications

* Coastal erosion monitoring
* Shoreline retreat assessment
* Climate change adaptation planning
* Sea-level rise vulnerability studies
* Coastal infrastructure risk analysis
* Protected area management

---

## Disclaimer

This project is an independent implementation inspired by the methodologies described in the USGS Digital Shoreline Analysis System (DSAS). It is not affiliated with, endorsed by, or maintained by the United States Geological Survey (USGS).

Users should independently validate results before applying them to regulatory or engineering decisions.

---

## Author

**Julieta Martin Tempestti**
Environmental Engineer | GIS & Coastal Change Analysis

---

## Citation

If you use this script in academic or professional work, please cite both:

* Himmelstoss et al. (2021), DSAS v5.1 User Guide.
* This repository.

