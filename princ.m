% Returns the principal value (between -pi and +pi) of an angle

function princ_val = princ(angle)

phi = mod(angle, 2*pi);
if phi < pi
    princ_val = phi;
else
    princ_val = phi - 2*pi;
end