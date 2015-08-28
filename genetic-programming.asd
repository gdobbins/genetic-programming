(in-package :asdf-user)

(defsystem "genetic-programming"
    :description "genetic-programming: genetic programming algorithms and other stochastic search methods."
    :version "0.0.1"
    :author "Graham Dobbins <gcdobbin@ncsu.edu>"
    :licence "GPLv3"
    :defsystem-depends-on ("alexandria" "lparallel")
    :components ((:file "gene")))
