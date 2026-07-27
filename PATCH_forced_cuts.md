# Patch "forced cuts" — tagli imposti dall'esterno (satellite)

Versione: `satellite_cuts` = `cuts_irri_test` (base + blackout irrigazione) **+** tagli forzati.

## 1. Perché

Nella versione standard il calendario dei tagli è deciso **offline** da `cropcoef`, una volta
sola, per combinazione *stazione meteo × coltura*, e memorizzato nelle tabelle `pheno/pheno_N/*.dat`.
Durante la simulazione IdrAgra si limita a consultarle:

```fortran
crop_pars_mat%k_cb(i,j) = info_pheno(ws_idx(i,j))%k_cb%tab(doy_s, soil_use%mat(i,j))
```

Il taglio, per il modello, **è** il crollo del Kcb dentro quella tabella. Di conseguenza:

* i tagli non possono variare da poligono a poligono (la tabella è condivisa da migliaia di celle);
* i tagli non possono reagire a nulla di ciò che accade nel bilancio idrico giornaliero.

Questa patch **disaccoppia la forma della curva dal calendario dei tagli**:

* `cropcoef` resta invariato e continua a produrre la *forma* della ricrescita (curva già calibrata);
* il *quando* tagliare diventa una proprietà **per cella**, decisa nel loop giornaliero.

## 2. Come funziona

Il puntatore fenologico `doy_s` viene **ri-ancorato**: quando una cella riceve un taglio imposto,
`doy_s` salta all'inizio della ricrescita di riferimento e poi avanza di un giorno al giorno.
La curva percorsa è quindi quella originale di `cropcoef`, solo fatta ripartire sulla data imposta.

* **Prima del primo taglio imposto dell'anno** il green-up primaverile resta quello vero, ma il
  puntatore viene **trattenuto sul plateau di maturità** (`regrow_start - 1`) e non può entrare
  nel taglio della curva di riferimento. Senza questo blocco il calendario GDD produrrebbe un
  taglio spurio prima di quello satellitare — e, se cadesse pochi giorni prima, lo
  **mascherebbe** del tutto (il Kcb sarebbe già crollato).
* **Dopo la fine della stagione di riferimento** (`season_end`, ultimo giorno di crescita
  attiva nella curva originale) l'ancoraggio viene **rilasciato**: la cella torna sulla
  curva originale, quindi senescenza autunnale e dormienza invernale restano quelle vere.
  Senza questo rilascio la coltura resterebbe verde fino a dicembre.
* **Se il taglio successivo tarda**, il puntatore si **ferma sul plateau di maturità**
  (`regrow_end`) invece di proseguire nella tabella e produrre un taglio non richiesto.
* La finestra di ricrescita (`regrow_start`, `regrow_end`) è ricavata automaticamente dalla
  curva di riferimento di ogni *stazione × coltura*: dal giorno di un taglio fino al giorno
  prima del taglio successivo.

## 3. Opt-in — nessuna regressione

Il meccanismo si attiva **solo** se esiste il file `geodata/forced_cuts.txt`.
Senza quel file il comportamento è identico a quello della versione precedente
(gli argomenti aggiunti sono `optional` in Fortran).

## 4. File di input

`{simulazione}/geodata/forced_cuts.txt` — scritto da IdragraTools in fase di export.

```
<numero di righe>
row	col	year	doy
45	42	2012	162
45	42	2012	207
...
```

* `row`, `col`: indici 1-based, stessa convenzione di `cells.txt` e dei raster `.asc`.
* `doy`: giorno dell'anno del taglio (le date satellitari hanno già l'offset di -7 giorni applicato
  a monte, nello script di inferenza `cuts_recognizer`).

Le celle **non** elencate nel file continuano a usare il normale calendario GDD.

## 5. File sorgente modificati

| File | Modifica |
|---|---|
| `mod_parameters.f90` | nuovo type `forced_cut` (row, col, year, doy) |
| `cli_watsources.f90` | nuova subroutine `open_forced_cuts` (gemella di `open_irrigation_blackout`) |
| `mod_crop_phenology.f90` | nuova subroutine `compute_regrow_window`; `populate_crop_pars_matrices` accetta 3 argomenti opzionali e ri-ancora `doy_s` |
| `cli_simulation_manager.f90` | lettura del file, calcolo della finestra di ricrescita a ogni anno, aggiornamento giornaliero di `fc_last_cut`, passaggio degli argomenti |
| `cli_read_parameter.f90` | **disattivati i 4 prompt interattivi** (`read *`), vedi §5.1 |

### 5.1 Prompt interattivi disattivati

Il codice originale si ferma in attesa di un INVIO in quattro punti puramente informativi
(righe 152, 273, 294, 321), il più fastidioso dei quali scatta quando la cartella `simout`
esiste già:

```
The directory simout\ already exists and will be updated
 <enter> to continue
```

Lanciato da QGIS o da uno script non c'è nessuno stdin che risponda, quindi il modello resta
bloccato a tempo indefinito senza dare errore. Le `read *` sono state commentate; i messaggi
informativi restano. Non è più necessario cancellare `simout` prima di ogni run.

## 6. Build

```
make cleanall && make
```

(oppure Ctrl+Shift+B in VS Code). Per un eseguibile di debug con backtrace: nel `Makefile`
decommentare la riga GFFLAGS di debug e commentare quella release, sempre preceduto da
`make cleanall`.

Con `-verbose` / `-v` il modello stampa il numero di record di `forced_cuts.txt` letti e,
in debug, l'elenco completo.

## 7. Verifica consigliata dopo la compilazione

1. Lanciare una simulazione **senza** `forced_cuts.txt` e confrontarla con l'output della
   versione precedente: devono coincidere (test di non-regressione).
2. Lanciare con il file su poche celle note (es. i 4 campi-stazione) e controllare in
   `simout/{anno}_cell_{row}_{col}.csv` che la colonna `kcb` crolli esattamente nei giorni
   indicati in `forced_cuts.txt`.

## 8. Sviluppo futuro reso possibile

Ora che la decisione del taglio è nel loop giornaliero, diventa semplice condizionarla anche
allo stress idrico reale (`hks`, calcolato in `calculate_water_stresses`) — l'intervento
"gate su hks" discusso per i tagli fantasma.
