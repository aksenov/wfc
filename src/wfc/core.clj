(ns wfc.core
  "Wave function collapse"
  (:require 
   [wfc.sample :as sample]
   [clojure.set :as set]
            [clojure.spec.alpha :as s]
            [clojure.math :as math]
            ))


(set! *warn-on-reflection* true)
(set! *unchecked-math* :warn-on-boxed)



;; sample loading
;; TODO in another module

;; sample analysis

(s/def :wfc/tile keyword?)
(s/def :wfc/sample-row (s/coll-of :wfc/tile :kind vector?))
(s/def :wfc/sample (s/coll-of :wfc/sample-row :kind vector?))


(defn entropy
  "Shannon's entropy"
  [variants weights]
  (let [tile-weights (map #(double (get weights %)) variants)
        tile-weight-sum (double (reduce + tile-weights))]
    (-
     (math/log tile-weight-sum)
     (/ (double (reduce + (map (fn [^double x] (* x (math/log x))) tile-weights)))
        tile-weight-sum))))

(defrecord Tile [r c collapsed? variants entropy index])

;; TODO get weights, variants from world and create new tile as world tile
(defn new-tile
  "Create new tile"
  [r c variants weights]
  (->Tile r c false variants (entropy variants weights) [r c])
  #_{:r r :c c :collapsed? false :variants variants :entropy (entropy variants weights)})

(defn collapse-tile
  "Collapse tile to certain variant"
  [tile variant]
  (if (:collapsed? tile)
    tile
    (assoc tile
           :variants variant
           :collapsed? true
           :entropy nil)))


(defn get-tile
  "Get tile from the `world` in position [`r` `c`]"
  ([world [r c]]
   (get-tile world r c))
  ([world r c]
   (get-in world [:state r c])))


(defn lowest-entropy-tile
  "Search for the tile with lowest entropy"
  [world]
  (if (:lowest-entropy world)
    (get-tile world (:lowest-entropy-coords world))
    (reduce
     (fn [res el]
       (if (or (nil? res) (< ^double (:entropy el) ^double (:entropy res)))
         el
         res))
     nil
     (filter
      :entropy
      (flatten
       (:state world))))))

;; TODO make position matter
(defn new-world
  "Build new `world`"
  [r c sample-analysis]
  (let [{:keys [adjacency weights superposition]} sample-analysis
        tile (new-tile 0 0 superposition weights)]
    {:state (into [] (for [r' (range r)]
                       (into [] (for [c' (range c)]
                                  (new-tile r' c' superposition weights))))
                  #_(repeat r (into [] (repeat c tile))))
     :rows r
     :cols c
     :collapsed? false
     :adjacency adjacency
     :weights weights
     :lowest-entropy (entropy superposition weights)
     :lowest-entropy-coords [(rand-int r) (rand-int c)]
     :superposition superposition}))


(defn adjacent-pos 
  "Get coordinates of adjacents tiles left, right, up and down."
  ([world r c]
   (adjacent-pos (:rows  world) (:cols world) r c))
  ([^long rows ^long cols ^long r ^long c]
   (let [left (dec c)
         right (inc c)
         up (dec r)
         down (inc r)]
     {:left (when (<= 0 left) [r left])
      :right (when (< right cols) [r right])
      :up (when (<= 0 up) [up c])
      :down (when (< down rows) [down c])})))


(defn adjacent-coords
  [^long r ^long c ^long rows ^long cols]
  (let [left (dec c)
        right (inc c)
        up (dec r)
        down (inc r)]
    [(when (<= 0 left) [r left])
     (when (< right cols) [r right])
     (when (<= 0 up) [up c])
     (when (< down rows) [down c])]))


(defn update-lowest-entropy
  "Remember lowest entropy and tile with it."
  ([world [r c] force?]
   (update-lowest-entropy world r c force?))
  ([world r c force?]
   (let [current-tile (get-tile world r c)
         ^double tile-entropy (:entropy current-tile)
         ^double lowest-entropy (:lowest-entropy world)]
     (cond
       force? 
       (assoc world
              :lowest-entropy tile-entropy 
              :lowest-entropy-coords [r c])
       
       (nil? tile-entropy) world
       
       (and (nil? lowest-entropy) (some? tile-entropy))
       (assoc world
              :lowest-entropy tile-entropy
              :lowest-entropy-coords [r c])
       
       (< tile-entropy lowest-entropy)
       (assoc world
              :lowest-entropy tile-entropy
              :lowest-entropy-coords [r c])
       :else world))))


(defn- rand-weighted-choice
  "Random key from weight map."
  [weights]
  (let [wsum (reduce + (vals weights))
        r (long (rand-int wsum))]
    (loop [w (seq weights)
           s 0]
      (let [[k ^long v] (first w)]
        (if (or (empty? w)
                (<= r (+ v s)))
          k
          (recur (rest w)
                 (+ v s)))))))


(defn random-weighted-variant
  "Select variant from possible `variants` according to weight."
  [variants weights]
  (rand-weighted-choice (select-keys weights variants)))

(defn zero-entropy?
  [tile]
  (when tile
    (when (set? (:variants tile))
      (= 1 (count (:variants tile))))
    #_(when-let [^double e (:entropy tile)]
        (<= e 0))))


(defn reduce-adjacent-variants
  "Reduce variants of adjacent tiles"
  [world adj-pos adj-tiles weights] 
  (reduce
   (fn [w dir]
     (if-let [dir-coords (get adj-pos dir)]
       (let [[dr dc] dir-coords
             tile (get-tile world dr dc)
             {:keys [variants collapsed?]} tile]
         (if (or collapsed? (= 1 (count variants)))
           w
           (let [dir-variants (get adj-tiles dir)
                 new-variants (set/intersection variants dir-variants)]
             (if (empty? new-variants) 
               (reduced nil)
               (-> w
                   (update-in [:state dr dc] assoc
                              :variants new-variants
                              :entropy (entropy new-variants weights))
                   (update-lowest-entropy dr dc false))))))
       w))
   world
   (shuffle [:left :right :up :down])))

(defn collapse
  "Collapse cell and reduce varability of it's adjacent tiles."
  [world r c variants]
  (loop [vs variants]
    (if (empty? vs)
      world ;; FIXME too permissive
      (let [variant (random-weighted-variant vs (:weights world))
            {:keys [weights adjacency]} world
            adj-pos (adjacent-pos world r c)
            adj-tiles (get adjacency variant)
            slightly-collapsed-world (-> world
                                         (update-in [:state r c] collapse-tile variant)
                                         (update-lowest-entropy r c true))
            world' (reduce-adjacent-variants slightly-collapsed-world adj-pos adj-tiles weights)
           ]
        (if world'
          world'
          (recur (set/difference vs #{variant})))))))


;; TODO progressive world weight

(defn collapse-adjacent
  ([world [r c] tile-variants]
   (collapse-adjacent world r c tile-variants))
  ([world r c tile-variants]
   (let [;v (if (set? variant) (first variant) variant)
         world'  (collapse world r c tile-variants)
         [left right up down] (adjacent-coords r c (:rows world') (:cols world'))]
     (reduce
      (fn [w dir-coords]
        (let [dir-tile (when dir-coords (get-tile w dir-coords))]
          (if (zero-entropy? dir-tile) 
            (collapse-adjacent w dir-coords (:variants dir-tile))
            w)))
      world'
      (shuffle [left right up down])))))


(defn collapse-world
  "Collapse whole world to the certanity"
  ([world]
   (collapse-world world (constantly true)))
  ([world continue-fn]
   (loop [world' world
          last-tile nil]
     (if (continue-fn world')
       (let [le-tile (lowest-entropy-tile world')]
         (if-not (or (nil? le-tile) (= le-tile last-tile))
           (recur (collapse-adjacent world' (:index le-tile) (:variants le-tile))
                  le-tile)
           (do (println le-tile)
               world')))
       world'))))


(defn print-collapsed
  [tile]
  (if (:uncertain? tile)
    (print (str (name :🟥)))
    (print (str (name (:variants tile))))))


(defn print-variants
  [tile]
  ;(print (str " " (count (:variants tile)) " "))
  (let [c (count (:variants tile))]
    (print (case c
             0 "0️⃣"
             1 "1️⃣"
             2 "2️⃣"
             3 "3️⃣"
             4 "4️⃣"
             5 "5️⃣"
             6 "6️⃣"
             "⬜"))))


(defn print-tile
  [tile]
  (if (:uncertain? tile)
    (print "❌")
    (if (:collapsed? tile)
      (print-collapsed tile)
      (print-variants tile))))


(defn print-world
  [world]
  (doseq [r (range (:rows world))]
    (doseq [c (range (:cols world))]
      (print-tile (get-tile world r c)))
    (println ":")))


(defn entropy
  "Shannon's entropy"
  [variants weights]
  (let [tile-weights (map #(get weights %) variants)
        ^double tile-weight-sum (reduce + tile-weights)]
    (-
     (Math/log tile-weight-sum)
     (/ (double (reduce + (map (fn [^double x] (* x (Math/log x))) tile-weights)))
        tile-weight-sum))))

(let [{:keys [weights]}
      (sample/analyze-sample [[:🟫 :🟫 :🟫 :🟫]
                       [:🟩 :🟫 :🟫 :🟩]
                       [:🟦 :🟩 :🟩 :🟦]
                       [:🟦 :🟦 :🟦 :🟦]])]
  (entropy [:🟦] weights))
(comment
  (time
   (let [s1 [[:🟩 :🌳 :🟩 :🟩 :🌳  :🟩 :🟩 :🟩]
             [:🟩 :🟩 :🟩 :🟩 :🟡  :🟡 :🟨 :🟩]
             [:🌳 :🌳 :🌳 :🟡 :🟦  :🟦 :🟨 :🌳]
             [:🌲 :🌲 :🌳 :🟡 :🟦  :🟦 :🏨 :🌳]
             [:🌳 :🌲 :🟩 :🟩 :🟡  :🟦 :🟨 :🌲]
             [:🌳 :🌳 :🟩 :🟡 :🟦  :🟦 :🟨 :🌳]
             [:🌾 :🌳 :🟩 :🟡 :🟦  :🟨 :🌾 :🟩]
             [:🟩 :🌳 :🏡 :🟡 :🟨  :🟨 :🟩 :🟩]
             [:🟩 :🟩 :🌳 :🟩 :🌳  :🟩 :🟩 :🟩]]
         s2 [[:🟫 :🟫 :🟫 :🟫]
             [:🟩 :🟫 :🟫 :🟩]
             [:🟦 :🟩 :🟩 :🟦]
             [:🟦 :🟦 :🟦 :🟦]]
         s3 [[:🌾 :🌾]
             [:🌾 :🌾]]
         s4 [[:🟩 :🟩 :🟩 :🟩]
             [:🟩 :🟦 :🟦 :🟦]
             [:🟩 :🟦 :🌴 :🟦]
             [:🟩 :🟦 :🟦 :🟦]]
         s6 [[:🌳 :🌳 :🌳 :🌳]
             [:🌳 :🌲 :🌲 :🌳]
             [:🌳 :🌲 :🌲 :🌳]
             [:🌳 :🌳 :🌳 :🌳]]
         s7 [[:🟦 :🟦 :🟦]
             [:🟦 :🏝 :🟦]
             [:🟦 :🟦 :🟦]]
         s10 [[:🟩 :🟩 :🟩 :🟩 :🟩 :🟩]
              [:🟩 :🟫 :🟫 :🟫 :🟫 :🟩]
              [:🟩 :🟫 :🟨 :🟨 :🟫 :🟩]
              [:🟩 :🟫 :🟨 :🟨 :🟫 :🟩]
              [:🟩 :🟫 :🟫 :🟫 :🟫 :🟩]
              [:🟩 :🟩 :🟩 :🟩 :🟩 :🟩]]
         s5 [[:A :A :A :A]
             [:A :A :A :A]
             [:A :A :A :A]
             [:A :C :C :A]
             [:C :B :B :C]
             [:C :B :B :C]
             [:C :B :B :C]]
         w (new-world 25 25 (sample/analyze-samples [s1 s6 s7]))]
     (-> w
         ;(collapse-adjacent 0 0 [:🟫])
         ;(collapse-adjacent 1 1 [:🟨])
         ;(collapse-adjacent 2 2 [:🟦])
       ;(update-lowest-entropy 3 3)
         (collapse-world #(do
                            ;(println (:lowest-entropy-coords %))
                            ;(print-world %)
                            %

                            true))
      ;(collapse-adjacent 0 3 :🏡)
      ;(lowest-entropy-tile)
         (print-world)
         (println))
  ;(collapse-adjacent w 1 1 :🟩)
  ;(lowest-entropy-tile )
     )))

;; TODO select with adjacency weights
;; TODO indexed search - X - use indexes inside tiles