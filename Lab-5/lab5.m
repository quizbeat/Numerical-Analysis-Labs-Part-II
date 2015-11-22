function [] = lab5()

fi_0 = @(x) (0);
fi_l = @(x) (0);
psi = @(x) (sin(2*pi.*x));

h = 0.01;
l_max = 2;
t_max = 2;
tao = 1;

X = 0:h:l_max;
T = 0:h:t_max;
            
u(1,:) = psi(X);
            
for k = 2:size(T)
    u(k,1) = fi_0(X(k));
    u(k,end) = fi_l(X(k));
    for j = 2:size(X)-1
        u(k,j) = u(k-1,j) - (tao/h^2) * (u(k-1,j+1) - 2*u(k-1,j) + u(k-1,j-1));
    end
end

plotmatrix(u)

end

