close all
clear all


data_dir = 'Z:\share\MNEs_from andrei\MNEs\';
data_files = dir([data_dir, 'st_750*.mat']);

color_range = [-.15, .15];

height = 16;
width = 20;

figure()
k=1; kmin=1; kmax=length(data_files); hk=loop_slider_n(k,kmin,kmax,1);
while true
    if ~ishandle(hk)
        break
    end
    k = round(get(hk, 'Value'));
    
    temp_handle = subplot(1,1,1);
    cla(temp_handle)

    load([data_dir, data_files(k).name])

    [V, D] = eig(reshape(J, 320, 320));

    [D, I] = sort(diag(D));
    V = V(:,I);
    Ndim = length(D);
    
    subplot(1, 2, 2)
    hold off
    plot(D, '.')
    hold on
    m = mean(D);
    s = std(D);
    std_factor = 2;
    plot([0,Ndim], [m-std_factor*s, m-std_factor*s], 'k--')
    plot([0,Ndim], [m+std_factor*s, m+std_factor*s], 'k--')
    title(['eigenvalues of ', data_files(k).name],'interpreter','none')
    
    num_features = sum(D>m+std_factor*s | D<m-std_factor*s);
    t = ceil(sqrt(num_features));
    t2 = ceil(num_features / t);
    
    counter = 1;
    for j = (1:sum(D>m+std_factor*s))-1
        subplot(t2, t*2, counter)
        imagesc(reshape(V(:,end-j), height, width), color_range)
        colormap(gray)
        counter = counter + 1;
        if mod(counter, t*2) / (t*2) > .5
            counter = counter + t;
        end
    end
    for j = 1:sum(D<m-std_factor*s)
        subplot(t2, t*2, counter)
        imagesc(reshape(V(:,j), height, width), color_range)
        colormap(gray)
        if j == 1
            title('-')
        end
        counter = counter + 1;
        if mod(counter, t*2) / (t*2) > .5
            counter = counter + t;
        end
    end
    uiwait;
end