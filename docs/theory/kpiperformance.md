# Key Performance indicators

Key Performance Indicators (KPIs) are essential metrics used to evaluate the performance of geothermal projects. In pyfastwell, KPIs are calculated based on the results of the DoubletCalc and economic model calculations. 
The KPIs provide insights into the efficiency, profitability, and sustainability of geothermal energy projects.

The KPI have been subdivided into the following categories:



## technical KPIs

Reservoir and well performance indicators are calculated based on the results of the DoubletCalc calculation

- **Reservoir Temperature [C]**: The temperature of the geothermal fluid at the reservoir, which is crucial for determining the energy potential.
- **Reservoir Depth [m]**: The depth below surface of the top of the reservoir
- **Thickness [m]**: The thickness of the reservoir
- **Net hydraulic Transmissivity [Dm]**: A measure of the reservoir's ability to transmit fluids, expressed in Darcymeters (Dm). 
It is calculated as the product of permeability and net thickness of the reservoir.
- **Production Temperature [C]**:  production temperature at the surface, 
typically 1-2 degrees C temperature drop compared to Reservoir temperature
- **Injection Temperature [C]**: injection temperature at the surface r.
- **Flow Rate [m3/h]**: The rate at which geothermal fluid is produced from the reservoir.
- **Pressure Drop [bar]**: The pressure  required to drive the thermal loop. The doublet model takes into account hydraulic resistance of the reservoir,
friction in tubing and thermosyphon effects.
- **Power[MW]**: Power production, which is the net power output of the geothermal plant after accounting for conversion efficiency and parasitic losses.
- **COP[-]**: The system Coefficient of Performance (COP),  which is the ratio of useful power provided to the energy consumed by the ESP and heat pump system.


## economic KPIs

- **LCOE [€cts/kWh]**: The Unit Technical Cost or socalled Levelized Cost of Energy, which is the average cost per unit of energy produced over the lifetime of the project, 
expressed in euro cts per kilowatt-hour.
- **NPV [million€]**: The net present value of the project, which is the difference between 
the present value of cash inflows and outflows over the project's lifetime. In pythermogis it is calculated as:

