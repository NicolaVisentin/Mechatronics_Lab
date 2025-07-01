function mouseClickCallback(~,~)

    % Definisci variabili globali per mantenere i punti tra le chiamate

    global x_points y_points

    % Se le variabili non sono ancora state inizializzate, inizializzale

    if isempty(x_points)
        x_points = [];
        y_points = [];
    end

    % Ottieni il punto dove è stato cliccato

    p = get(gca, 'Currentpoint');  % Ottieni la posizione del click rispetto agli assi
    x = p(1,1);  % La posizione x
    y = p(1,2);  % La posizione y

    % Aggiungi il punto all'elenco

    x_points = [x_points, x];
    y_points = [y_points, y];

    % Traccia i punti

    plot(x_points, y_points, 'ro-', 'MarkerFaceColor','r');

    % Verifica se è stato cliccato con il tasto destro

    if strcmp(get(gcf, 'SelectionType'), 'alt')  % 'alt' è il tasto destro

        % Salva la traiettoria nel workspace
        assignin('base', 'x_points_acquired', x_points);
        assignin('base', 'y_points_acquired', y_points);

        % Disabilita la callback per fermare l'esecuzione
        set(gcf, 'WindowButtonDownFcn', []);
        uiresume;

    end

end