classdef NALab5
    
    methods(Static)
        
        function [output] = solve(l,T,fi_0,fi_l,psi)
            h = 0.001;
            
            u_prev = psi(0:h:l);
            
            for t = 1:h:T
               for x = 1:h:l
                  u(row,1) = fi_0(x,t);
                  u(row,end) = fi_l(x,t);
               end
            end
        end
        
        function [] = part1()
            % input data
            fi_0 = @(x) (0);
            fi_l = @(x) (0);
            psi = @(x) (sin(2*pi.*x));
            
            out = solve(l,T,fi_0,fi_l,psi);
            
        end
        
    end
    
end

