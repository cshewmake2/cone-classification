function G = gauss2d(x, y, mu, sigma)
    exponent = ((x-mu(1)).^2 + (y-mu(2)).^2)./(2*sigma^2);
    amplitude = 1 / (sigma * sqrt(2*pi));  
    G = amplitude * exp(-exponent);
end