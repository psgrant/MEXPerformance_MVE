function phi_list = eval_pade(z, k, d)
[s, z] = scaled_arg(z);
Id = speye(size(z));
[phi_vals{1:k}] = scaled_pade(z, Id, d, k);     % (d,d)-PadÈ approx

if s > 0,
    [phi_list{1:s+1}] = square_pade(z, Id, s, phi_vals);
else
    phi_list = {phi_vals};
end

end