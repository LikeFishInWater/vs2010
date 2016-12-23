function int_var = gf_transform_2_int(gf_var, order, primpoly)
for int_var = 0:2^order-1
    if gf(int_var, order, primpoly) == gf_var
        break;
    end
end
end

