function setFigureMarkerSize(size)

kids = get(gca, 'Children');

for (kk = kids)
    try
        set(kk, 'MarkerSize', size);
    end
end