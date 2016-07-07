<?xml version="1.0"?>
<ogs6>
<coupling>
    <P algorithm="Serial" convergence="FemFunctionConvergenceCheck" max_itr="1" epsilon="1e-4">
        <out>GAS_PRESSURE</out>
        <out>CAPILLARY_PRESSURE</out>
		<out>TEMPERATURE</out>
        <problems>
            <M name="CMP_NonIso_CapPressureForm" type="CMP_NonIso_CapPressureForm">
                <out>GAS_PRESSURE</out>
				<out>CAPILLARY_PRESSURE</out>
				<out>TEMPERATURE</out>
            </M>
        </problems>
    </P>
</coupling>
</ogs6> 

