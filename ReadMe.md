# Elite Präsentation
Am Institut für Theoretische Chemie der Universität Stuttgart.

## Interaktives Online Notebook
In diesem Online Notebook könnt ihr vorgefertigte Skripte und Visualisierungen testen.
Klick einfach auf die Badge.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/KoehnLab/ElitePresentation/HEAD?labpath=tutorial.ipynb)

## Zum selbst Installieren
Um das Notebook auf eurem Laptop oder PC auszuprobieren, müsst ihr zuerst das Projekt herunterladen und ein paar Softwarepakete installieren.

### Download
Klickt weiter oben auf den grünen **Code** Button und wählt anschließend **Download Zip**.
Das Zip kann nach dem Herunterladen in einem Ordner eurer Wahl entpackt werden.
<details>
  <summary>Alternative (klick hier)</summary>
  Falls ihr euch mit Git auskennt, könnt ihr natürlich den <code>git clone</code> Befehl verwenden um den Sourcecode herunterzuladen.
</details>


### Python
#### Windows
Falls ihr es noch nicht installiert habt: [Python](https://www.python.org/downloads/)

#### Linux
Falls noch kein Python auf eurem System vorhanden ist, könnt ihr es mit Hilfe eures Packagemanagers auf einem Terminal installieren.
Da der Syntax leider von dem installierten OS abhängt zeigen wir hier nur eine (Ubuntu) der vielen Möglichkeiten.
```
sudo apt install python3 python3-pip
```

### Python Pakete
Mit der Installation von Python wird auch der Pythoneigene Packagemanager *pip* mitgeliefert.
Mit diesem können nun alle notwendigen Pakete installiert werden.
Alle nötigen Pakete sind in der Datei [requirements.txt](requirements.txt) aufgelistet.
In dem Terminal könnt ihr ein Paket zum Beispiel so
```
pip install numpy
```
installieren.
Um alle Pakete auf einmal zu installieren könnt ihr den Befehl
```
pip install -r requirements.txt
```
verwenden.

### Jupyter Lab
Um mit Jupyter Notebooks arbeiten zu können empfehlen wir [Jupyter Lab](https://jupyter.org/install) zu installieren.
```
pip install jupyterlab
```
Anschließend wollen wir *Jupyter Lab* in dem Projektordner starten.

#### Windows
Verwendet den *Explorer*, aufrufbar durch die Tastenkombination `Win+E` und navigiert zu dem entpackten oder geklonten Projekt.
Dort angekommen könnt ihr in der Pfadleiste (dort wo der momentane Ordnerpfad angezeigt wird) einfach den Befehl `powershell` eingeben und mit `Enter` bestätigen. Nun sollte sich ein Terminal öffnen.
Anschließend könnt ihr in das Terminal `jupyter-lab` eingeben und wieder mit `Enter` bestätigen.
Diese Eingabe sollte euren Browser öffnen und das Projekt in der Seitenleiste links anzeigen.

#### Linux
Verwendet ein Terminal um zu dem entpackten oder geklonten Projekt zu navigieren.
```
cd <Pfad_zum_projekt>
```
Anschließend führt den Befehl
```
jupyter-lab
```
aus und euer Browser sollte sich mit dem Projekt darin öffnen.

<details>
<summary>Alternative zu Jupyter Lab (klick hier)</summary>
Alternativ kann das Projekt auch direkt in einem Texteditor aufgerufen, ausgeführt und bearbeitet werden.
Hierfür könnt ihr zum Beispiel <a href="https://code.visualstudio.com/download">VScode</a> verwenden.
</details>
