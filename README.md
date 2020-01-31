# slice-img-ASTRALIS

analyze data collected by GCaMP calcium-sensor imaging in astrocytes in mouse brain slices.

1. Requires user to input basic information regarding:
   - mouse, 
   - slice, 
   - imaging conditions,
   - and all relevant experimental variables.
   
2. After full fov analysis using ImageJ macros, loads the data into Matlab for analysis.

3. After detection of calcium activity in video recordings using Matlab/ImageJ-based application (AQuA: https://github.com/yu-lab-vt/AQuA) for non-ROI based detection of calcium events, data is transferred to Matlab.

4. Analysis includes> 
   - event properties: number, size, duration, delay, propagation;
   - differences in events across experimental conditions;
   - pooling of data across subjects/slices.
