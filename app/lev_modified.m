function [OUT] = lev_modified(IN, m_view, current_operation_view, delay)
    
    
    len = length(IN);
    rho0 = get_R(IN, 0);
    rho = zeros(1,len-1);
    delta = zeros(1,len-1);
    nabla = zeros(1,len-1);
    a = zeros(len-1,len-1);
    b = zeros(len-1,len-1);
    OUT = zeros(len,len);

    for m = 1:len-1
        m_view.Text = sprintf("%d из %d",m,len-1);
        pause(delay);
        current_operation_view.Text = "Вычисляется дельта";
        for k = 1:m
            pause(delay);
            if m-k == 0
                delta(m) = delta(m) + get_R(IN, k) * 1;
            else
                delta(m) = delta(m) + get_R(IN, k) * a(m-1,m-k);
            end
        end

        current_operation_view.Text = "Вычисляется набла";
        for k = 1:m
            pause(delay);
            if m-k == 0
            nabla(m) = nabla(m) + get_R(IN, -k) * 1;
            else
                nabla(m) = nabla(m) + get_R(IN, -k) * b(m-1,m-k);
            end
        end
        
        current_operation_view.Text = "Вычисляются a,b,p";
        if m-1 ~= 0
            a(m,m) = -delta(m)/rho(m-1);
            b(m,m) = -nabla(m)/rho(m-1);
            rho(m) = rho(m-1) - (delta(m)*nabla(m)/rho(m-1));
        else
            a(m,m) = -delta(m)/rho0;
            b(m,m) = -nabla(m)/rho0;
            rho(m) = rho0 - (delta(m)*nabla(m)/rho0);
        end

        for k=1:m-1
            pause(delay);
            if m-1 ~= 0
                a(m,k) = a(m-1,k)+a(m,m)*b(m-1,m-k);
                b(m,k) = b(m-1,k)+b(m,m)*a(m-1,m-k);
         
            end

        end
        
        

       if m == len-1

        current_operation_view.Text = "Вычисляются значения обратной матрицы";
        OUT(1,1) = 1/rho(m);
        OUT(len,len) = 1/rho(m);

        for k=1:len-1
            pause(delay);
            OUT(1,k+1) = a(len-1,k)/rho(m);
            OUT(len-k, len) = a(len-1,k)/rho(m);
        end
        
        for k=1:len-1
            pause(delay);
            OUT(k+1,1) = b(len-1,k)/rho(m);
            OUT(len, len - k) = b(len-1,k)/rho(m);
        end


            for j=1:len-2
                for k=1:len-2
                  pause(delay);
                  OUT(j+1,k+1) = OUT(j,k) + ((1/rho(m))*(a(m,j)*b(m,k)-a(m,m-k+1)*b(m,m-j+1)));
                end
            end
       end

       OUT(len,len) = OUT(1,1);
       
    end
    
end