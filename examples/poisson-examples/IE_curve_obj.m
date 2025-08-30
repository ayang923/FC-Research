classdef IE_curve_obj < handle
    %IE_CURVE_OBJ Curve obj for integral equation solver
    
    properties
        l_1
        l_2
        l_1_prime
        l_2_prime
        l_1_dprime
        l_2_dprime
        
        n
        
        w
        w_prime
        
        s_mesh
        theta_mesh
        
        curve_idx
        next_curve
    end
    
    methods
        function obj = IE_curve_obj(curve, curve_idx, cfac, p)
            %IE_CURVE_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj.l_1 = curve.l_1;
            obj.l_2 = curve.l_2;
            obj.l_1_prime = curve.l_1_prime;
            obj.l_2_prime = curve.l_2_prime;
            obj.l_1_dprime = curve.l_1_dprime;
            obj.l_2_dprime = curve.l_2_dprime;
            
            obj.n = ceil((curve.n-1)*cfac);
            
            obj.curve_idx = curve_idx;
            curve_idx
            
            v = @(s) (1/p-1/2)*(1-2*s).^3+1/p*(2*s-1)+1/2;
            v_prime = @(s) 2/p-6*(1/p-1/2)*(1-2*s).^2;

            obj.w = @(s) v(s).^p./(v(s).^p+v(1-s).^p);
            obj.w_prime = @(s) p*((v_prime(s).*v(s).^(p-1))./(v(s).^p+v(1-s).^p)-(v(s).^(p-1).*v_prime(s)-v(1-s).^(p-1).*v_prime(1-s)).*v(s).^p./(v(s).^p+v(1-s).^p).^2);
            
            ds = 1/obj.n;
            obj.s_mesh = linspace(0, 1-ds, obj.n);
            
            obj.theta_mesh = obj.w(obj.s_mesh);
        end
        
    end
end


