classdef IE_curve_seq_obj < handle
    %IE_CURVE_SEQ_OBJ Obj for integral equation solver
    
    properties
        n_curves
        first_curve
        last_curve
    end
    
    methods
        function obj = IE_curve_seq_obj(curve_seq, cfac, p)
            %IE_CURVE_SEQ_OBJ Construct an instance of this class
            
            obj.n_curves = curve_seq.n_curves;
            obj.first_curve = IE_curve_obj(curve_seq.first_curve, 1, cfac, p);
            
            curr_curve = curve_seq.first_curve;
            curr_ie_curve = obj.first_curve;
            for i = 2:obj.n_curves
                curr_curve= curr_curve.next_curve;
                curr_ie_curve.next_curve = IE_curve_obj(curr_curve, i, cfac, p);
                curr_ie_curve = curr_ie_curve.next_curve;
            end
            obj.last_curve = curr_ie_curve;
            curr_ie_curve.next_curve = obj.first_curve;
        end
        
    end
end

