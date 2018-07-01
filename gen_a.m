function a = gen_a(M,Delta,theta)

for f = 1:length(Delta)
    for k = 1:M
        a(k,f) = exp(1i*(2*pi)*((k-1)*Delta(f))*sind(theta));
    end
end

end