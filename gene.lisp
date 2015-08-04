(in-package :cl-user)
;(eval-when (:compile-toplevel :load-toplevel :execute)
  ;(ql:quickload :series)
 ; (ql:quickload :alexandria)
  ;(ql:quickload :lparallel))
(defpackage :genetic-programming
  (:use :common-lisp :alexandria :lparallel)
  (:export :gene-prog :make-population :population-genes :population-gene-fitness-list :population-gene-adjusted-fitness-list :population-number-of-genes :population-average-fitness :population-maximum-fitness :population-minimum-fitness :population-best-gene-so-far :population-best-gene-so-far-fitness :make-functions :functions-function-list :functions-arg-num-list :simulated-annealing :genetic-algorithm :make-program :program-lisp-function :program-genome :program-fitness-history :program-value-history :program-running-fitness :program-running-value :program-most-recent-output :program-number-of-parents :program-number-of-parents-left :program-number-of-parents-right :program-creation-method :program-most-recent-fitness :program-most-recent-value :simple-running-fitness-tabulator :make-program-pool :program-pool-programs :program-pool-size :program-pool-functions :program-rank :program-pool-terminals :program-pool-program-evaluator :program-pool-fitness-function :program-pool-value-function :program-pool-running-fitness-tabulator :program-pool-running-value-tabulator :create-random-full-gene-program :create-simple-child :create-simple-comparison-child :roullete-wheel-chooser :fill-pool-with-simple-children :fill-pool-proportionally-with-simple-children :fill-pool-with-simple-comparison-children :fill-pool-proportionally-with-simple-comparison-children :fill-pool-rank-proportionally-with-simple-comparison-children :feed-program-unknown-data :feed-pool-unknown-data :tabulate-pool-running-fitness :cull-bottom-of-pool :feed-program-known-data :feed-pool-known-data :fill-program-pool-randomly :make-mod-2-program-evaluator :straight-running-fitness-tabulator :print-program-pool-information :simple-program-pool-crossover :delete-program-pool-fitness-history :fill-pool-with-random-children :fill-pool-with-top-of-other-pool :create-simple-sub-pools :fill-sub-pools-randomly :feed-sub-pools-unknown-data :cross-sub-pools :collapse-simple-sub-pools :cull-bottom-of-pool-list :print-program-pool-list-information :fill-sub-pools-proportionally-with-simple-children))
(in-package :genetic-programming)

;(defclass gene ()
;  ((gene :accessor gene :initarg :gene)
;   (raw-fitness :accessor raw-fitness)
;   (adjusted-fitness :accessor adjusted-fitness)))

(defvar *number-of-threads* 4)
(setf lparallel:*kernel* (lparallel:make-kernel *number-of-threads*))



(defstruct population
  (genes nil :type (or cons null))
  (gene-fitness-list nil :type (or cons null))
  (gene-adjusted-fitness-list nil :type (or cons null))
  (number-of-genes 0 :type (and fixnum (real 0)))
  (average-fitness 0 :type number)
  (maximum-fitness 0 :type number)
  (minimum-fitness 0 :type number)
  (best-gene-so-far nil :type (or cons null))
  (best-gene-so-far-fitness 0 :type number))

(defstruct functions
  (function-list nil :type (or cons null))
  (arg-num-list nil :type (or cons null))
  (base-gene-copy-number 1 :type (and fixnum (real 1))))

;(defmacro concurrent-lists-do (lists &body body)
;  (let* ((syms (loop for i in lists collect (gensym))) (len (length syms)))
;    `(do ,(loop for i from 0 to (1- len) collect `(,(nth i syms) ,(nth i lists) (rest ,(nth i syms))))
;	 ((or ,@(loop for i from 0 to (1- len) collect `(null ,(nth i syms)))))
;       (let ,(loop for i from 0 to (1- len) collect `(,(nth i lists) ,(nth i syms)))
;       ,@body))))

(defun tree-replace (tree find max-depth min-depth replace &rest vars)
  (declare (type cons tree) (type fixnum max-depth min-depth) (type (or null cons) vars) (type function replace))
  (loop for i in tree collect
       (if (and (consp i) (<= 0 max-depth) (not #1=(find-if (lambda (x) (equalp i x)) (last vars))))
	   (apply 'tree-replace i find (1- max-depth) (1- min-depth) replace vars)
	   (if (or (equal i find) (and (<= 0 min-depth) #1#) (and (>= 0 max-depth) (consp i) (not #1#)))
	       (apply replace vars)
	       i))))

(defun tree-find-not (tree find)
  (declare (type cons tree))
  (let ((rep t))
    (declare (type boolean rep))
    (loop for i in tree do
	 (if rep
	     (if (consp i)
		 (setf rep (tree-find-not i find))
		 (if (equal i find)
		     (setf rep nil)))))
    rep))

(defun tree-depth (tree &optional (depth 0) (max-depth 0) (ret nil))
  (declare (type (or cons null) tree ret) (type (and fixnum (real 0)) depth max-depth))
  (if tree
      (if (consp (car tree))
	  (tree-depth (cdar tree) (1+ depth) (max max-depth (1+ depth)) (cons (list depth (cdr tree)) ret))
	  (tree-depth (cdr tree) depth max-depth ret))
      (if ret
	  (tree-depth (second (first ret)) (first (first ret)) max-depth (rest ret))
	  max-depth)))

(defun tree-atom-size (tree &optional (functions nil) (atoms 0) (ret nil))
  "Counts the number of atoms in a tree. Does not count () expressions as atoms, skips the function slot of the list, i.e. the first position of every list if functions is nil. Subtract one if gene is not wrapped in a list."
  (declare (type (or cons null) tree ret) (type boolean functions) (type (and fixnum (real 0)) atoms))
  (if tree
      (if (consp (car tree))
	  (tree-atom-size (cdar tree) functions (if functions (1+ atoms) atoms) (cons (cdr tree) ret))
	  (tree-atom-size (cdr tree) functions (1+ atoms) ret))
      (if ret
	  (tree-atom-size (first ret) functions atoms (rest ret))
	  atoms)))

(defun node-list (gene &optional (acc nil) (ret nil))
  (declare (type (or cons null) gene acc ret))
  (if gene
      (if (consp (car gene))
	  (node-list (cdar gene) (cons gene acc) (cons (cdr gene) ret))
	  (node-list (cdr gene) (cons gene acc) ret))
      (if ret
	  (node-list (first ret) acc (rest ret))
	  acc)))

(defun create-random-gene (functions terminals &optional (max-depth 6) (min-depth 2))
  (declare (type functions functions) (type cons terminals))
  (declare (type (and fixnum (real -1)) max-depth min-depth))
  (let* ((len-fun (length #1=(functions-function-list functions))) (len-term (length terminals))
	 (len-fun-mult (if (> len-term len-fun) (round (/ len-term len-fun)) 1))
	 (new-len-fun (* len-fun len-fun-mult))
	 (len-term-mult (if (> len-fun len-term) (round (/ len-fun len-term)) 1))
	 (len (+ (* len-fun len-fun-mult) (* len-term len-term-mult))) gene)
    (declare (type (and fixnum (real 0)) len-fun len-term len-term-mult len-fun-mult new-len-fun len) (type (or null cons) gene))
    (flet ((create-new-sexp (num functions terminals)
	     (declare (type fixnum num) (type functions functions) (type cons terminals))
		 (if (< num new-len-fun)
		     (append (list (nth (mod num len-fun) #1#)) (make-list (nth (mod num len-fun) (functions-arg-num-list functions)) :initial-element :sexp))
		     (nth (mod (- num new-len-fun) len-term) terminals))))
      (setf gene (create-new-sexp (random len-fun) functions terminals))
      (do () ((tree-find-not gene :sexp) gene)
	(setf gene (tree-replace gene :sexp max-depth min-depth (lambda (x y z) (declare (type (and fixnum (real 0)) x)) (create-new-sexp (random x) y z)) len functions terminals))))))

(defun create-random-full-gene (functions terminals &optional (max-depth 6) (min-depth 2))
  (declare (type functions functions) (type cons terminals) (type fixnum max-depth min-depth))
  (loop for i from 1 to (functions-base-gene-copy-number functions) collect
       (create-random-gene functions terminals max-depth min-depth)))

(defun initialize-population (current-population functions terminals generation)
  (declare (type population current-population) (type cons terminals) (type functions functions) (type fixnum generation))
  (assert (= 0 generation))
  (if (< (length #1=(population-genes current-population)) #2=(population-number-of-genes current-population))
      (tagbody (setf (population-genes current-population) (append (list (create-random-full-gene functions terminals)) (population-genes current-population))) (initialize-population current-population functions terminals generation))
;      (initialize-population (make-population :genes (append (list (create-random-gene functions terminals)) #1#) :number-of-genes #2#) functions terminals generation)
      current-population))

(defun pre-test-population (current-population functions terminals pt-function)
  (declare (type population current-population) (type cons terminals) (type functions functions) (type function pt-function))
  (dolist (i (population-genes current-population))
    (funcall pt-function (nth 0 i) 0 (the boolean t))
    (dotimes (j (functions-base-gene-copy-number functions))
      (do () ((the boolean (funcall pt-function (nth j i) j (the boolean nil))))
	(let ((temp (create-random-gene functions terminals)))
	  (declare (type cons temp))
	  (rplaca (nth j i) (car temp))
	  (rplacd (nth j i) (cdr temp)))))))

(defun adjust-to-linear-fitness (fitness-list max-fit min-fit avg-fit fit-list-length &optional (max-fit-multiplier 2))
  (declare (type cons fitness-list) (type real max-fit min-fit avg-fit) (type (and fixnum (real 0)) fit-list-length max-fit-multiplier))
  (if (= max-fit avg-fit) (make-list fit-list-length :initial-element 1.0)
      (let* ((a (/ (- max-fit-multiplier 1) (- max-fit avg-fit))) (b (* -1 (1+ (* a avg-fit)))) (result #1=(pmapcar (lambda (x) (declare (type number x)) (+ b (* a x))) :size fit-list-length fitness-list)))
	(declare (type number a b) (type cons result))
;	(setf a (/ (- max-fit-multiplier 1) (- max-fit avg-fit)) b (* -1 (1+ (* a avg-fit))))
;	#1=(setf result (pmapcar (lambda (x) (declare (type number x)) (+ b (* a x))) fitness-list)
	(if (psome #'minusp :size fit-list-length result)
	    (if (= avg-fit min-fit)
		(setf result (make-list fit-list-length :initial-element 1.0))
		(setf a (/ 1 (- avg-fit min-fit)) b (* -1 a min-fit) result #1#)))
	(pmap-into result (lambda (x) (declare (type real x)) (coerce x 'single-float)) :size fit-list-length result))))

(defun update-population-fitness (current-population fitness-function)
  (declare (type function fitness-function) (type population current-population))
  (labels ((upf-best-gene-helper (gene-list fitness-list best-gene-fitness)
	     (declare (type cons gene-list fitness-list) (type real best-gene-fitness))
	     (if (= (the real (car fitness-list)) best-gene-fitness)
		 (values (car gene-list) (car fitness-list))
		 (upf-best-gene-helper (cdr gene-list) (cdr fitness-list) best-gene-fitness))))
    (pmap-into #1=(population-gene-fitness-list current-population) fitness-function :size (population-number-of-genes current-population) (population-genes current-population))
    (setf (population-average-fitness current-population) (/ (preduce #'+ #1#) (population-number-of-genes current-population)))
    (setf (population-minimum-fitness current-population) (preduce #'min #1#))
    (setf (population-maximum-fitness current-population) (preduce #'max #1#))
    (setf (population-gene-adjusted-fitness-list current-population) (adjust-to-linear-fitness (population-gene-fitness-list current-population) (population-maximum-fitness current-population) (population-minimum-fitness current-population) (population-average-fitness current-population) (population-number-of-genes current-population)))
    (if (> (population-maximum-fitness current-population) (population-best-gene-so-far-fitness current-population))
	(multiple-value-bind (gene fitness) (upf-best-gene-helper (population-genes current-population) (population-gene-fitness-list current-population) (population-maximum-fitness current-population))
	  (setf (population-best-gene-so-far current-population) gene (population-best-gene-so-far-fitness current-population) fitness)))))

(defun cmp-helper (gene-list fitness-list selector-number)
  (declare (type single-float selector-number) (type cons gene-list fitness-list) (optimize (speed 3) (safety 0) (compilation-speed 0)))
  (if (> (the single-float (car fitness-list)) selector-number)
      (car gene-list)
      (if (cdr gene-list)
	  (cmp-helper (cdr gene-list) (cdr fitness-list) (- selector-number (the single-float (car fitness-list))))
	  (cmp-helper gene-list fitness-list (1- selector-number)))))

(defun create-mating-pool (current-population)
  (declare (type population current-population))
  (plet ((total-fitness (reduce (lambda (x y) (declare (type single-float x y) (optimize (speed 3) (safety 0))) (handler-case (+ x y) (floating-point-overflow () most-positive-single-float))) (population-gene-adjusted-fitness-list current-population))))
    (declare (type single-float total-fitness))
    (make-population :number-of-genes (population-number-of-genes current-population) :best-gene-so-far (population-best-gene-so-far current-population) :best-gene-so-far-fitness (population-best-gene-so-far-fitness current-population)
		     :average-fitness (population-average-fitness current-population)
		     :maximum-fitness (population-maximum-fitness current-population)
		     :minimum-fitness (population-minimum-fitness current-population)
		     :gene-fitness-list (population-gene-fitness-list current-population)
		     :gene-adjusted-fitness-list (population-gene-adjusted-fitness-list current-population)
		     :genes (pmapcar (lambda (x) (declare (ignore x)) (copy-tree (cmp-helper (population-genes current-population) (population-gene-adjusted-fitness-list current-population) (random total-fitness)))) (population-genes current-population)))))

;		     (loop for i from 1 to (population-number-of-genes current-population) collect (copy-tree (cmp-helper (population-genes current-population) (population-gene-adjusted-fitness-list current-population) (random total-fitness)))))))
(defun bmp-helper (gene-list functions terminals pairs-remaining mutate-percentage &optional (max-depth 50))
  (declare (type (and fixnum (real 0)) pairs-remaining max-depth) (type cons terminals) (type (or null cons) gene-list) (type (single-float 0.0 1.0) mutate-percentage)
	   (type functions functions))
  (if (<= 2 pairs-remaining)
      (let ((base-gene-number (random (functions-base-gene-copy-number functions))))
	(plet ((sub-one (random-elt (node-list (list (nth base-gene-number (first gene-list)))))) (sub-two (random-elt (node-list (list (nth base-gene-number (second gene-list)))))))
	  (declare (type cons sub-one sub-two))
	  (rotatef (car sub-one) (car sub-two))
	  (plet ((tree-depth-of-second-gene-list (tree-depth (second gene-list))))
	    (if (< max-depth (the fixnum (tree-depth (first gene-list))))
		(if #1=(< max-depth (the fixnum tree-depth-of-second-gene-list))
		    (rotatef (car sub-one) (car sub-two))
		    (setf (car sub-one) (car sub-two)))
		(if #1#
		    (setf (car sub-two) (car sub-one)))))
	  (bmp-helper (cddr gene-list) functions terminals (- pairs-remaining 2) mutate-percentage max-depth)))
      (if gene-list
	  (let ((n (length gene-list)))
	    (declare (type fixnum n))
	    (dotimes (i (the fixnum (round (+ (* n mutate-percentage (- 1 mutate-percentage) (gaussian-random)) (* n mutate-percentage)))))
	      (plet ((sub-gene (random-elt gene-list)))
		(plet ((sub-one (random-elt (node-list sub-gene))) (sub-two (list (create-random-gene functions terminals -1 1))))
		  (rotatef (car sub-one) (car sub-two))
		  (if (< max-depth (the fixnum (tree-depth sub-gene)))
		      (rotatef (car sub-one) (car sub-two))))))))))

(defun breed-mating-pool (mating-pool functions terminals &optional (breed-percentage 1.0) (mutate-percentage 0.0))
  (declare (type population mating-pool) (type functions functions) (type cons terminals) (type (single-float 0.0 1.0) breed-percentage mutate-percentage))
  (bmp-helper (population-genes mating-pool) functions terminals (round (* (population-number-of-genes mating-pool) breed-percentage)) mutate-percentage)
    mating-pool)

(declaim (inline update-history-report))
(defun update-history-report (current-population generation print-freq)
  (declare (type fixnum generation print-freq))
  (if (zerop (mod generation print-freq))
      (format t "~&--------------------------------------------------~%Generation: ~d     Stardate: ~{~d-~d-~d  ~d:~2,'0d:~2,'0d~}~%Maximum Fitness: ~d   Minimum Fitness: ~d   Average Fitness: ~d   Best Fitness: ~d   "
	      generation (nreverse (subseq (multiple-value-list (get-decoded-time)) 0 6)) (population-maximum-fitness current-population) (population-minimum-fitness current-population) (population-average-fitness current-population) (population-best-gene-so-far-fitness current-population))))

(let ((current-population (make-population)) (mating-pool (make-population))
      (terminals nil) (functions (make-functions)) (print-freq 1)
      (fitness-function #'(lambda (x) (eval (car x))))
      (generation 0))
  (declare (type fixnum generation print-freq) (type (or null cons) terminals) (type function fitness-function) (type functions functions) (type population current-population mating-pool))
  (defun gene-prog (&key (operation :next-generation) (generations 50);for progressing to the next evolution
		      (delete-duplicates nil) (population-size 100) (init-functions functions) (init-terminals terminals) (fit-function fitness-function) (print-frequency print-freq);for during initialization, whether to check for duplicates
		      (pt-function nil)
		      (file nil)) ;for saving and loading
    (declare (type (and fixnum (real 0)) generations population-size print-frequency) (type keyword operation) (type boolean delete-duplicates) (type functions init-functions) (type cons init-terminals) (type function fit-function) (type (or null (simple-array character (*)))))
    (case operation
      (:next-generation
       (update-population-fitness current-population fit-function)
       (setf mating-pool (create-mating-pool current-population))
       (setf current-population (breed-mating-pool mating-pool functions terminals))
       (update-history-report current-population generation print-freq)
       (setf generation (the (and fixnum (real 0)) (1+ generation))))
      (:best-of-generation
       (labels ((bog-helper (gene-list fitness-list max-fitness)
		  (if (= max-fitness (car fitness-list))
		      (car gene-list)
		      (bog-helper (cdr gene-list) (cdr fitness-list) max-fitness))))
	 (update-population-fitness current-population fit-function)
	 (values (bog-helper (population-genes current-population) (population-gene-fitness-list current-population) (population-maximum-fitness current-population)) (population-maximum-fitness current-population))))
      (:best-of-generation-continue
       (setf mating-pool (create-mating-pool current-population))
       (setf current-population (breed-mating-pool mating-pool functions terminals))
       (update-history-report current-population generation print-freq)
       (setf generation (the (and fixnum (real 0)) (1+ generation))))
      (:next-generations
       (dotimes (i generations)
	 (update-population-fitness current-population fitness-function)
	 (setf mating-pool (create-mating-pool current-population))
	 (setf current-population (breed-mating-pool mating-pool functions terminals))
	 (update-history-report current-population generation print-freq)
	 (setf generation (the (and fixnum (real 0)) (1+ generation)))))
       (:initialize
       (setf terminals init-terminals functions init-functions
	     (population-number-of-genes current-population) population-size
	     (population-number-of-genes mating-pool) population-size)
       (setf fitness-function fit-function)
       (setf print-freq print-frequency)
       (initialize-population current-population functions terminals generation)
       (setf (population-gene-fitness-list current-population) (make-list (population-number-of-genes current-population) :initial-element 0))
       (if delete-duplicates (funcall 'gene-prog :operation :delete-initial-duplicates
				      :population-size population-size
				      :init-functions init-functions
				      :init-terminals init-terminals
				      :delete-duplicates t)))
       (:pre-test
	(pre-test-population current-population functions terminals pt-function))
      (:enlarge-population
       (setf (population-number-of-genes current-population) population-size
	     (population-number-of-genes mating-pool) population-size)
       (initialize-population current-population functions terminals 0)
       (setf (population-gene-fitness-list current-population) (make-list (population-number-of-genes current-population) :initial-element 0)))
      (:delete-initial-duplicates
       (let ((old-population-size (length (population-genes current-population))))
       (setf #1=(population-genes current-population) (delete-duplicates #1# :test #'tree-equal))
       (funcall 'gene-prog :operation :initialize
		:population-size old-population-size
		:init-functions init-functions
		:init-terminals init-terminals
		:print-frequency print-frequency
		:delete-duplicates (if (= old-population-size (length #1#)) nil t))))
      (:save (with-open-file (out file :direction :output :if-exists :supersede)
	       (prin1 current-population out)
	       (prin1 mating-pool out)
	       (prin1 terminals out)
	       (prin1 functions out)
	       (prin1 (list generation) out)
	       (prin1 (list print-freq) out)))
      (:load (with-open-file (in file :direction :input)
	       (setf current-population (read in))
	       (setf mating-pool (read in))
	       (setf terminals (read in))
	       (setf functions (read in))
	       (setf generation (car (read in)))
	       (setf print-freq (car (read in)))))
      (:change-terminals (setf terminals init-terminals))
      (:set-fit-function (setf fitness-function fit-function))
      (:mating-pool mating-pool)
      (:fitness-function fitness-function)
      (:terminals terminals)
      (:functions functions)
      (:generation generation)
      (:print-freq print-freq)
      (:current-population current-population))))

(defun boltzmann (e-i e-i+1 temp &optional (k 1.0d0))
  (declare (type real e-i e-i+1) (type (double-float 0.0d0) temp k))
  (exp (the double-float (/ (abs (- e-i e-i+1)) -1 k temp))))

(defun simulated-annealing (initial-value initial-temperature energy-function step-function &optional (min-t 2.0d-6) (mu-T 1.003d0) (k 1.0d0))
  (declare (type function energy-function step-function) (type double-float initial-temperature min-t mu-t k))
  (let ((value initial-value) (value-energy (funcall energy-function initial-value)))
    (do* ((candidate-value (funcall step-function value) (funcall step-function value)) (candidate-value-energy (funcall energy-function candidate-value) (funcall energy-function candidate-value)) (temp initial-temperature (/ temp mu-T))) ((< temp min-t) (if (< candidate-value-energy value-energy) (values candidate-value candidate-value-energy) (values value value-energy)))
      (declare (type double-float temp))
      (if (or (< candidate-value-energy value-energy) (< (random 1.0d0) (the double-float (boltzmann value-energy candidate-value-energy temp k))))
	  (setf value candidate-value value-energy candidate-value-energy)))))

(export 'genetic-algorithm)
(defun genetic-algorithm (terminal-list gene-length population-size generations fitness-function &key (mutate 0.01) (verbose nil) (result :default) (greedy nil))
  (declare (type cons terminal-list) (type function fitness-function) (type (and fixnum (real 0)) gene-length population-size generations) (type keyword result) (type boolean verbose greedy))
  (labels ((gene-maker (terminal-list gene-length &optional (acc nil) (len (length acc)))
	     (declare (type cons terminal-list) (type (or null cons) acc) (type (and fixnum (real 0)) gene-length len))
	     (if (= len gene-length)
		 acc
		 (gene-maker terminal-list gene-length (cons (random-elt terminal-list) acc) (1+ len))))
	   (gene-list-maker (terminal-list gene-length population-size &optional (acc nil) (len (length acc)))
	     (declare (type cons terminal-list) (type (and fixnum (real 0)) gene-length population-size len) (type (or null cons) acc))
	     (if (= len population-size)
		 acc
		 (gene-list-maker terminal-list gene-length population-size (cons (gene-maker terminal-list gene-length nil 0) acc) (1+ len))))
	   (cmg-helper (current-genes fitness-list fitness-residual)
	     (declare (type cons current-genes fitness-list) (type real fitness-residual))
	     (if (< fitness-residual (the real (first fitness-list)))
		 (copy-list (first current-genes))
		 (if (rest current-genes)
		     (cmg-helper (rest current-genes) (rest fitness-list) (- fitness-residual (the real (first fitness-list))))
		     (cmg-helper current-genes fitness-list (- fitness-residual 1)))))
	   (create-mating-genes (current-genes fitness-list)
	     (declare (type cons current-genes fitness-list))
	     (if greedy
		 (append #4=(list (copy-list (nth (position (apply #'max fitness-list) fitness-list) current-genes)))
			 #4# #4#
			 (let ((total-fitness (apply #'+ fitness-list)))
			   (loop for i from 1 to (- population-size 4) collect (cmg-helper current-genes fitness-list (random total-fitness)))) #4#)
		 (let ((total-fitness (apply #'+ fitness-list)))
			   (loop for i from 1 to population-size collect (cmg-helper current-genes fitness-list (random total-fitness))))))
	   (breed-mating-genes (mating-genes &optional (fin mating-genes))
	     (declare (type (or cons null) mating-genes fin))
	     (if (second mating-genes)
		 (let* ((place (random (1- gene-length))) (sub-one (nthcdr place (first mating-genes))) (sub-two (nthcdr place (second mating-genes))))
		 (rotatef (cdr sub-one) (cdr sub-two))
		 (breed-mating-genes (rest (rest mating-genes)) fin))
		 fin)))
    (let ((current-genes (gene-list-maker terminal-list gene-length population-size))
	  (current-genes-fitness (make-list population-size :initial-element 0))
	  (current-genes-adjusted-fitness (make-list population-size :initial-element 0.0))
	  (mating-genes (make-list population-size)))
      (declare (type cons current-genes current-genes-fitness current-genes-adjusted-fitness mating-genes))
      (dotimes (i generations)
	#1=(pmap-into current-genes-fitness fitness-function :size population-size current-genes)
	(setf current-genes-adjusted-fitness (adjust-to-linear-fitness current-genes-fitness
								      (apply #'max current-genes-fitness)
								      (apply #'min current-genes-fitness)
								      (/ (apply #'+ current-genes-fitness) population-size)
								      population-size 2))
	(setf mating-genes (create-mating-genes current-genes current-genes-adjusted-fitness))
	(if verbose
	    (format t "~&Generation: ~d   Best of Generation: ~a   Fitness: ~d" i
		    (nth (position #2=(apply #'max current-genes-fitness) current-genes-fitness) current-genes)
		    #2#))
	(setf current-genes (breed-mating-genes mating-genes))
	(pmap-into current-genes (lambda (x) (declare (type (or null cons) x)) (if (zerop (random (round (* mutate population-size))))
						(let ((place (random gene-length)))
						  (append (subseq x 0 place) (list (random-elt terminal-list)) (subseq x (1+ place))))
						x)) current-genes))
      (case result
	(:default #1# (values (nth (position (apply #'max current-genes-fitness) current-genes-fitness) current-genes) #2#))))))
	
;;;;;EUGENE

(defun default-program-lisp-function (&rest args)
  (error "default-program-lisp-function was called with args: ~a" args))


(defstruct program
  (lisp-function #'default-program-lisp-function :type function)
  (genome nil :type (or cons null))
  (fitness-history nil :type (or null cons))
;  (value-history nil :type (or null cons))
  (running-fitness 0 :type real)
;  (running-value 0 :type real)
;  (most-recent-fitness 0 :type real)
  (most-recent-output nil :type (or null cons))
  (rank 0 :type (and fixnum (real 0)))
;  (genome-history-left nil :type (or null cons))
;  (genome-history-right nil :type (or null cons))
  (number-of-parents-left 0 :type (integer 0))
  (number-of-parents-right 0 :type (integer 0))
  (creation-method :default :type keyword :read-only t))

(declaim (inline program-number-of-parents))
(defun program-number-of-parents (program-instance)
  (declare (type program program-instance))
  (the (integer 0) (+ (program-number-of-parents-left program-instance) (program-number-of-parents-right program-instance))))

(declaim (inline program-most-recent-fitness))
(defun program-most-recent-fitness (program-instance)
  (declare (type program program-instance))
  (the real (first (program-fitness-history program-instance))))

;(declaim (inline program-most-recent-value))
;(defun program-most-recent-value (program-instance)
;  (declare (type program program-instance))
;  (the real (first (program-value-history program-instance))))

(defun straight-running-fitness-tabulator (fitness-history)
  (declare (type cons fitness-history))
  (loop for i of-type real in fitness-history for j of-type (and fixnum (real 1)) from 1 summing i into sum of-type real finally (return (/ sum j))))

(defun simple-running-fitness-tabulator (fitness-history)
  (declare (type cons fitness-history))
  (loop for i of-type real in fitness-history for j of-type (and fixnum (real 1)) from 1 to 100 summing (/ i j) into sum of-type real summing (/ 1 j) into recip of-type (and rational (real 1)) finally (return (/ sum recip))))

(defstruct program-pool
  (programs nil :type (or null cons))
  (size 0 :type (and fixnum (real 0)))
  (functions nil :type functions :read-only t)
  (terminals nil :type (or null cons) :read-only t)
  (program-evaluator nil :type function :read-only t)
  (fitness-function #1=(lambda (x y) (declare (ignore y) (type cons x)) (car x)) :type function)
;  (value-function #1# :type function)
  (running-fitness-tabulator #'simple-running-fitness-tabulator :type function))
;  (running-value-tabulator #'straight-running-fitness-tabulator :type function))

(defun delete-program-pool-fitness-history (program-pool)
  (declare (type program-pool program-pool))
  (loop for prgm of-type (or null program) in (program-pool-programs program-pool)
       if prgm do
       (setf (program-fitness-history prgm) nil); (program-value-history prgm) nil)
     end
     finally (return program-pool)))

(defmacro make-mod-2-program-evaluator (variable-list)
  (declare (optimize (speed 3) (space 0) (compilation-speed 0) (safety 0)))
  `(eval `(lambda ,,variable-list (declare (ignorable ,@,variable-list) (optimize speed compilation-speed)) (mapcar (lambda (z) (mod z 2)) (list ,@gene)))))

(defun create-random-full-gene-program (functions terminals program-evaluator &optional (max-depth 6) (min-depth 2))
  (declare (type functions functions) (type cons terminals) (type fixnum max-depth min-depth) (type function program-evaluator))
  (let ((temp
	 (make-program :genome (create-random-full-gene functions terminals max-depth min-depth)
		       :creation-method :crfgp)))
    (declare (type program temp))
    (setf (program-lisp-function temp) (funcall program-evaluator (program-genome temp)))
    temp))

(defun create-simple-child (parent-1 parent-2 program-evaluator &optional (base-genome-length (min (length (program-genome parent-1)) (length (program-genome parent-2)))) (creation-method :csc))
  (declare (type program parent-1 parent-2) (type function program-evaluator) (type keyword creation-method) (type (and fixnum (real 0)) base-genome-length))
  (let ((temp (make-program :genome (copy-tree (program-genome parent-1))
;			    :genome-history-left (cons (program-genome parent-1) (program-genome-history-left parent-1))
;			    :genome-history-right (cons (program-genome parent-2) (program-genome-history parent-2))
			    :number-of-parents-left (the (integer 0) (+ 1 (program-number-of-parents parent-1)))
			    :number-of-parents-right (the (integer 0) (+ 1 (program-number-of-parents parent-2)))
			    :creation-method creation-method)))
    (declare (type program temp))
    (let ((base-place (random base-genome-length)))
      (declare (type (and fixnum (real 0)) base-place))
      (setf (car (random-elt (node-list (list (nth base-place (program-genome temp)))))) (car (random-elt (node-list (list (nth base-place (program-genome parent-2))))))))
    (setf (program-lisp-function temp) (funcall program-evaluator (program-genome temp)))
    temp))

(defun create-simple-comparison-child (parent-1 parent-2 program-evaluator ideal-output &optional (base-genome-length (min (length (program-genome parent-1)) (length (program-genome parent-2)))) (creation-method :cscc))
  (declare (type program parent-1 parent-2) (type function program-evaluator) (type keyword creation-method) (type (and fixnum (real 0)) base-genome-length) (type (or null cons) ideal-output))
  (if ideal-output (assert (= (length ideal-output) base-genome-length)))
  (let ((temp (make-program :genome (copy-tree (program-genome parent-1))
;			    :genome-history-left (cons (program-genome parent-1) (program-genome-history parent-1) (program-genome-history parent-2) (list (program-genome parent-2)))
			    :number-of-parents-left (the (integer 0) (+ 1 (program-number-of-parents parent-1)))
			    :number-of-parents-right (the (integer 0) (+ 1 (program-number-of-parents parent-2)))
			    :creation-method creation-method)))
    (declare (type program temp))
    (let* ((parent-1-places (mapcar #'equalp (program-most-recent-output parent-1) ideal-output))
	   (possibility-places (loop for i of-type boolean in (mapcar (lambda (x y) (declare (type boolean x y)) (if x nil (if y t))) parent-1-places (mapcar #'equalp (program-most-recent-output parent-2) ideal-output)) for j of-type (and fixnum (real 0)) from 0 nconc (if i (list j))))
	   (base-place (if possibility-places (random-elt possibility-places) (if (some (lambda (x) (declare (type boolean x)) x) parent-1-places) (random-elt (loop for i of-type boolean in parent-1-places for j of-type (and fixnum (real 0)) from 0 nconc (if i (list j)))) (random base-genome-length)))))
      (declare (type (or null cons) possibility-places parent-1-places) (type (and fixnum (real 0)) base-place))
      (setf (car (random-elt (node-list (list (nth base-place (program-genome temp)))))) (car (random-elt (node-list (list (nth base-place (program-genome parent-2))))))))
    (setf (program-lisp-function temp) (funcall program-evaluator (program-genome temp)))
    temp))

(defmacro linear-fill-pool-with-children (program-pool child-maker creation-method &optional extra-input)
  `(do ((children (do ((lst (program-pool-programs ,program-pool) (rest lst)) (count 0 (1+ count)))
		      ((null (car lst)) (if (not (assert (< (program-pool-size ,program-pool) (* count 2)))) lst))
		    (declare (type (and fixnum (real 0)) count))) (rest children))
	(parents (program-pool-programs ,program-pool) (rest parents)))
       ((null children) program-pool)
     ,(if extra-input
	  `(rplaca children (,child-maker (first parents) (second parents) (program-pool-program-evaluator ,program-pool) ,extra-input (functions-base-gene-copy-number (program-pool-functions ,program-pool)) ,creation-method))
	  `(rplaca children (,child-maker (first parents) (second parents) (program-pool-program-evaluator ,program-pool) (functions-base-gene-copy-number (program-pool-functions ,program-pool)) ,creation-method)))))

(defun roullete-wheel-chooser (choice-list weight-list &optional (total-weight (apply #'+ weight-list)))
  (declare (type cons choice-list weight-list) (type single-float total-weight))
  (labels ((rwc-helper (c-list w-list selector)
	     (declare (type cons c-list w-list) (type single-float selector))
	     (if (< selector (the single-float (car w-list)))
		 (car c-list)
		 (if (cdr w-list)
		     (rwc-helper (cdr c-list) (cdr w-list) (- selector (the single-float (car w-list))))))))
    (do ((temp #1=(rwc-helper choice-list weight-list (random total-weight)) #1#)) (temp temp))))

(defmacro proportional-fill-pool-with-children (program-pool child-maker creation-method &optional extra-input by-rank)
  (declare (type boolean by-rank))
  `(let* ((run-fit-list (loop for i of-type (or null program) in (program-pool-programs ,program-pool) nconc (if i (list ,(if by-rank '(program-rank i) '(program-running-fitness i))))))
	  (adj-fit-list (apply #'adjust-to-linear-fitness run-fit-list (loop for i of-type real in run-fit-list for j of-type (and fixnum (real 0)) from 1 summing i into sum of-type real maximizing i into max of-type real minimizing i into min of-type real finally (return (list max min (/ sum j) j 2)))))
	  (tot-fit (loop for i of-type single-float in adj-fit-list summing i into sum of-type single-float finally (return sum))))
     (do ((children (do ((lst (program-pool-programs ,program-pool) (rest lst)) (count 0 (1+ count)))
			((null (car lst)) lst)
		      (declare (type (and fixnum (real 0)) count))) (rest children)))
	 ((null children) program-pool)
       ,(if extra-input
	    `(rplaca children (,child-maker #1=(roullete-wheel-chooser (program-pool-programs ,program-pool) adj-fit-list tot-fit) #1# (program-pool-program-evaluator ,program-pool) ,extra-input (functions-base-gene-copy-number (program-pool-functions ,program-pool)) ,creation-method))
	    `(rplaca children (,child-maker #1# #1# (program-pool-program-evaluator ,program-pool) (functions-base-gene-copy-number (program-pool-functions ,program-pool)) ,creation-method))))))

(defun fill-pool-with-simple-children (program-pool)
  (declare (type program-pool program-pool))
  (do ((children (do ((lst (program-pool-programs program-pool) (rest lst)) (count 0 (1+ count))) ((null (car lst)) (if (not (assert (< (program-pool-size program-pool) (* count 2)))) lst)) (declare (type (and fixnum (real 0)) count))) (rest children))
       (parents (program-pool-programs program-pool) (rest parents)))
      ((null children) program-pool)
    (rplaca children (create-simple-child (first parents) (second parents) (program-pool-program-evaluator program-pool) (functions-base-gene-copy-number (program-pool-functions program-pool)) :csc-fpwsc))))

(defun fill-pool-proportionally-with-simple-children (program-pool)
  (declare (type program-pool program-pool))
  (proportional-fill-pool-with-children program-pool create-simple-child :csc-fppwsc))

(defun fill-pool-with-simple-comparison-children (program-pool ideal-output)
  (declare (type program-pool program-pool) (type (or null cons) ideal-output))
  (linear-fill-pool-with-children program-pool create-simple-comparison-child :cscc-fpwscc ideal-output))

(defun fill-pool-proportionally-with-simple-comparison-children (program-pool ideal-output)
  (declare (type program-pool program-pool) (type (or null cons) ideal-output))
  (proportional-fill-pool-with-children program-pool create-simple-comparison-child :cscc-fppwscc ideal-output))

(defun fill-pool-rank-proportionally-with-simple-comparison-children (program-pool ideal-output)
  (declare (type program-pool program-pool) (type (or null cons) ideal-output))
  (proportional-fill-pool-with-children program-pool create-simple-comparison-child :cscc-fprpwscc ideal-output t))

(defun fill-pool-with-random-children (program-pool &optional (max-depth 6) (min-depth 2))
  (declare (type program-pool program-pool) (type fixnum max-depth min-depth))
  (loop for prgms on (program-pool-programs program-pool)
       if (null (car prgms)) do
       (rplaca prgms
	       (create-random-full-gene-program (program-pool-functions program-pool)
						(program-pool-terminals program-pool)
						(program-pool-program-evaluator program-pool)
						max-depth min-depth)))
  program-pool)

(defun fill-pool-with-top-of-other-pool (program-pool other-pool)
  (declare (type program-pool program-pool other-pool))
  (loop for prgms on (program-pool-programs program-pool)
     with other-prgms of-type list = (program-pool-programs other-pool)
     if (null (car prgms)) do
       (rplaca prgms (car other-prgms)) and do
       (setf other-prgms (rest other-prgms)))
  program-pool)

(defun feed-program-unknown-data (program data fitness-function)
  (declare (type program program) (type (or null cons) data) (type function fitness-function))
  (let ((temp-output (apply (program-lisp-function program) data)))
    (setf #1=(program-fitness-history program) (nconc (list (funcall fitness-function temp-output data)) #1#))
    program))

(defun feed-pool-unknown-data (program-pool data)
  (declare (type program-pool program-pool) (type (or null cons) data))
  (dolist (i (program-pool-programs program-pool) program-pool)
    (feed-program-unknown-data i data (program-pool-fitness-function program-pool))))

(defun tabulate-pool-running-fitness (program-pool)
  (declare (type program-pool program-pool))
  (dolist (i (program-pool-programs program-pool) program-pool)
    (declare (type program i))
    (setf (program-running-fitness i) (funcall (program-pool-running-fitness-tabulator program-pool) (program-fitness-history i)))))

(defun cull-bottom-of-pool (program-pool &optional (cull-fraction 0.2))
  (declare (type program-pool program-pool) (type (single-float 0.0 1.0) cull-fraction))
  (setf #1=(program-pool-programs program-pool) (loop for prgm in (remove-duplicates (sort (program-pool-programs (tabulate-pool-running-fitness program-pool)) #'> :key #'program-running-fitness) :test #'tree-equal :key #'program-genome) for q from 1 to (round (* (- 1 cull-fraction) (program-pool-size program-pool))) collect prgm))
  (setf #1# (nconc #1# (make-list (the (and fixnum (real 0)) (- (program-pool-size program-pool) (length #1#))) :initial-element nil)))
  program-pool)

(defun feed-program-known-data (program data fitness-function ideal-output)
  (declare (type program program) (type (or null cons) data ideal-output) (type function fitness-function))
  (let ((temp-output (apply (program-lisp-function program) data)))
    (declare (type (or null cons) temp-output))
    (let ((temp-fitness (funcall fitness-function temp-output data)))
      (declare (type real temp-fitness))
      (setf #1=(program-fitness-history program) (nconc (list temp-fitness) #1#)
	    (program-most-recent-output program) temp-output
	    (program-rank program) (loop for i in temp-output for j in ideal-output sum (if (equalp i j) 1 0) into tally of-type (and fixnum (real 0)) finally (return tally))))
    program))

(defun feed-pool-known-data (program-pool data &optional ideal-output (cull-fraction 0.2) (child-creation-function #'fill-pool-proportionally-with-simple-comparison-children))
  (declare (type program-pool program-pool) (type (or null cons) data ideal-output) (type single-float cull-fraction) (type function child-creation-function))
  (funcall child-creation-function (dolist (i (program-pool-programs (cull-bottom-of-pool program-pool cull-fraction)) program-pool)
    (if i (feed-program-known-data i data (program-pool-fitness-function program-pool) ideal-output) (return program-pool))) ideal-output))

(defmacro fill-program-pool (program-pool fill-function)
  `(setf (program-pool-programs ,program-pool)
	 (loop for i from 1 to (program-pool-size ,program-pool) collect
	      (funcall ,fill-function))))

(defun fill-program-pool-randomly (program-pool &optional (max-depth 6) (min-depth 2))
  (declare (type program-pool program-pool) (type fixnum max-depth min-depth))
  (fill-program-pool program-pool (lambda () (create-random-full-gene-program (program-pool-functions program-pool)
								   (program-pool-terminals program-pool)
								   (program-pool-program-evaluator program-pool)
								   max-depth min-depth)))
  program-pool)

(defun print-program-pool-information (program-pool &optional (counter 0) extra-information)
  (declare (type program-pool program-pool) (type (integer 0) counter))
  (let ((first-program (first (program-pool-programs program-pool))))
    (declare (type program first-program))
    (format t "~&Time:                               ~{~4d-~2,'0d-~2,'0d  ~d:~2,'0d:~2,'0d~}~%Counter:                               ~d~%~a~%--------------------------------------------------------------------------------~%~%First program running fitness:                           ~10,3f~%~%First program rank:                                       ~10d~%First program number of parents:                 ~10d~%~%--------------------------------------------------------------------------------"
	    (nreverse (subseq (multiple-value-list (get-decoded-time)) 0 6))
	    counter
	    extra-information
	    (program-running-fitness first-program)
	    (program-rank first-program)
	    (program-number-of-parents first-program))
    program-pool))

(defun simple-program-pool-crossover (program-pool-1 program-pool-2)
  (declare (type program-pool program-pool-1 program-pool-2))
;  (assert (= (length (program-pool-programs program-pool-1)) (length (program-pool-programs program-pool-2))))
  (do ((p1 (program-pool-programs program-pool-1) (rest (rest p1))) (p2 (program-pool-programs program-pool-2) (rest (rest p2)))) ((null p1) (values program-pool-1 program-pool-2))
    (rotatef (car p1) (car p2))))

(defun create-simple-sub-pools (program-pool &optional sub-size)
  (declare (type program-pool program-pool))
  (let ((constants (remove-if-not #'constantp (program-pool-terminals program-pool)))
	(terms (remove-if #'constantp (program-pool-terminals program-pool))))
    (unless sub-size (setf sub-size (round (/ (* 4 (program-pool-size program-pool)) (length terms)))))
    (loop for term on terms
       collect (make-program-pool
		:size sub-size
		:functions (program-pool-functions program-pool)
		:terminals (append (list (car term) (if (cdr term) (cadr term) (car terms))) constants)
		:program-evaluator (program-pool-program-evaluator program-pool)
		:fitness-function (program-pool-fitness-function program-pool)
		:running-fitness-tabulator (program-pool-running-fitness-tabulator program-pool)))))

(defun collapse-simple-sub-pools (pool-list &optional size)
  (declare (type list pool-list) (type (or null number) size))
  (let ((in-size (if size size (/ (* (program-pool-size (first pool-list)) (length pool-list)) 4)))
	(len (length pool-list)))
    (make-program-pool
     :size in-size
     :programs (if (zerop (mod in-size len))
		   (apply #'append (mapcar (lambda (prg) (subseq (program-pool-programs prg) 0 (/ in-size len))) pool-list))
		   (apply #'append (mapcar (let ((cnt (mod in-size len))) (lambda (prg) (subseq (program-pool-programs prg) 0 (if (<= 0 (decf cnt)) (1+ (floor (/ in-size len))) (floor (/ in-size len)))))) pool-list)))
     :functions (program-pool-functions (first pool-list))
     :terminals (reduce #'union (mapcar #'program-pool-terminals pool-list))
     :program-evaluator (program-pool-program-evaluator (first pool-list))
     :fitness-function (program-pool-fitness-function (first pool-list))
     :running-fitness-tabulator (program-pool-running-fitness-tabulator (first pool-list)))))

(defmacro call-on-pool-list (func pool-list)
  `(mapcar ,func ,pool-list))

(defun fill-sub-pools-randomly (pool-list &optional (max-depth 6 max-depth-p) (min-depth 2 min-depth-p))
  (declare (type list pool-list))
  (if (or max-depth-p min-depth-p)
      (call-on-pool-list (lambda (pool) (fill-program-pool-randomly pool max-depth min-depth)) pool-list)
      (call-on-pool-list #'fill-program-pool-randomly pool-list)))

(defun feed-sub-pools-unknown-data (pool-list data)
  (declare (type list pool-list))
  (call-on-pool-list (lambda (pool) (feed-pool-unknown-data pool data)) pool-list))

(defun cull-bottom-of-pool-list (pool-list &optional (cull-fraction 0.2))
  (declare (type list pool-list))
  (call-on-pool-list (lambda (prg) (cull-bottom-of-pool prg cull-fraction)) pool-list))

(defun print-program-pool-list-information (pool-list &optional (counter 0))
  (declare (type list pool-list) (type number counter))
  (print-program-pool-information (first pool-list) counter)
  pool-list)

(defun fill-sub-pools-proportionally-with-simple-children (pool-list)
  (declare (type list pool-list))
  (call-on-pool-list #'fill-pool-proportionally-with-simple-children pool-list))

(defun cross-sub-pools (pool-list)
  (declare (type list pool-list))
  (rotate
   (loop for p1 in pool-list by #'cddr
      for p2 in (rest pool-list) by #'cddr
      nconcing (multiple-value-list (simple-program-pool-crossover p1 p2)) into temp of-type list
      finally (return (if (length= pool-list temp) temp (append temp (last pool-list)))))))
