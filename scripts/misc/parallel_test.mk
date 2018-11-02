# Teste.

all: ponga penha parra
	echo "$@ done!"


ponga: dagoberto
	./scripts/sleeper.sh 3 $@ "Foi3"

penha:
	./scripts/sleeper.sh 3 $@ "Foi3"

parra:
	./scripts/sleeper.sh 3 $@ "Foi3"


dagoberto:
	./scripts/sleeper.sh 5 $@ "Foi5"


print:
	echo "Only print."
